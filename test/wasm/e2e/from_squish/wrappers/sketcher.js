/**
 * Source-parity Playwright implementation of the Squish `Sketcher` wrapper.
 *
 * Method names deliberately follow
 * qa-ims/squish_modules/maestro/sketcher.py.  This wrapper is incremental:
 * a method is added only when its standalone WASM interaction has been
 * implemented with real browser input.
 */
import {
  clipboardText,
  clickMenuAction,
  clickPopupRow,
  clickWidget,
  drawingAreaCenter,
  exportText,
  focusCanvas,
  hoverPopupRow,
  isEmpty,
  loadStructureForTest,
  mouseClick,
  mouseDrag,
  openSketcher,
} from './sketcher_wasm.js';

const BUTTON_NAMES = {
  clear: 'clear_btn',
  clear_selection: 'clear_selection_btn',
  invert_selection: 'invert_selection_btn',
  move_rotate: 'move_rotate_btn',
  select_all: 'select_all_btn',
  undo: 'undo_btn',
};

const TOOL_NAMES = {
  C: 'c_btn',
  Cl: 'cl_btn',
  N: 'n_btn',
  O: 'o_btn',
  P: 'p_btn',
  S: 's_btn',
  F: 'f_btn',
  H: 'h_btn',
  move_rotate: 'move_rotate_btn',
  rect_btn: 'select_tool_btn',
};

const BOND_TOOL_NAMES = {
  1: 'single_bond_btn',
  2: 'bond_order_btn',
};

function v3000Block(text, name) {
  const begin = `M  V30 BEGIN ${name}`;
  const end = `M  V30 END ${name}`;
  const lines = String(text).split(/\r?\n/);
  const start = lines.indexOf(begin);
  const finish = lines.indexOf(end);
  if (start === -1 || finish === -1 || finish <= start) {
    throw new Error(`V3000 ${name} block not found`);
  }
  return lines.slice(start + 1, finish);
}

/** Parse the atom/bond information used by Squish get_structure_information(). */
function parseV3000Structure(text) {
  const atoms = v3000Block(text, 'ATOM').map((line) => {
    const fields = line.trim().split(/\s+/);
    const charge = Number((fields.find((field) => field.startsWith('CHG=')) || 'CHG=0').slice(4));
    return {
      index: Number(fields[2]),
      element: fields[3],
      x: Number(fields[4]),
      y: Number(fields[5]),
      charge,
    };
  });
  const bonds = v3000Block(text, 'BOND').map((line) => {
    const fields = line.trim().split(/\s+/);
    return {
      index: Number(fields[2]),
      order: Number(fields[3]),
      atom1: Number(fields[4]),
      atom2: Number(fields[5]),
      orientation: fields.find((field) => field.startsWith('CFG=')) || null,
    };
  });
  return { atoms, bonds };
}

/** Match Squish convert_to_sketcher_coordinates() for visible replay input. */
function sketcherCoordinates(structure) {
  const visibleAtoms = structure.atoms.filter(
    (atom) => !['A', 'AH', 'Q', 'QH', 'M', 'MH', 'X', 'XH'].includes(atom.element),
  );
  const center = visibleAtoms.reduce((sum, atom) => ({ x: sum.x + atom.x, y: sum.y + atom.y }), {
    x: 0,
    y: 0,
  });
  center.x /= visibleAtoms.length;
  center.y /= visibleAtoms.length;
  const byIndex = new Map(structure.atoms.map((atom) => [atom.index, atom]));
  const lengths = structure.bonds.map((bond) => {
    const first = byIndex.get(bond.atom1);
    const second = byIndex.get(bond.atom2);
    return Math.round(Math.hypot(second.x - first.x, second.y - first.y) * 1000) / 1000;
  });
  const counts = new Map();
  for (const length of lengths) counts.set(length, (counts.get(length) || 0) + 1);
  const commonLength = [...counts.entries()].sort((a, b) => b[1] - a[1])[0]?.[0] || 1.418;
  return {
    ...structure,
    atoms: structure.atoms.map((atom) => ({
      ...atom,
      x: ((atom.x - center.x) * 1.418 * 35) / commonLength,
      y: ((atom.y - center.y) * 1.418 * 35) / commonLength,
    })),
  };
}

const MORE_ACTION_NAMES = {
  add_explicit_hydrogens: 'Add Explicit Hydrogens',
  clear_selection: 'Clear Selection',
  copy_all_as: 'Copy All As',
  copy_all: 'Copy All',
  cut: 'Cut',
  fit_to_screen: 'Fit to Screen',
  flip_horizontal: 'Flip Horizontal',
  flip_vertical: 'Flip Vertical',
  invert_selection: 'Invert Selection',
  paste: 'Paste',
  redo: 'Redo',
  remove_explicit_hydrogens: 'Remove Explicit Hydrogens',
  select_all: 'Select All',
  undo: 'Undo',
};

const COPY_ALL_AS_NAMES = {
  cxsmi: 'Extended SMILES',
  inchikey: 'InChIKey',
  inchi: 'InChI',
  pdb: 'PDB',
  sdf: 'MDL SD V3000',
  smi: 'SMILES',
};

/** Browser equivalent of the source Squish `Sketcher` class. */
export class Sketcher {
  /** @param {import('@playwright/test').Page} page */
  constructor(page) {
    this.page = page;
    this.rebuild_structures = process.env.PLAYWRIGHT_REBUILD_STRUCTURES === '1';
    this.replay_tool = null;
    this.current_tool = null;
    this.replay_atom_indices = new Map();
    this.replay_bond_indices = new Map();
  }

  async open() {
    await openSketcher(this.page);
  }

  /**
   * Fixture setup only; never use this method in `tst_import_menu`.
   * Set PLAYWRIGHT_REBUILD_STRUCTURES=1 to replay Squish's visible
   * get_structure_information -> clear -> build_structure workflow.
   */
  async load_structure_for_test(text, { rebuild = this.rebuild_structures, sort = true } = {}) {
    await loadStructureForTest(this.page, text);
    if (!rebuild) return;
    const structure = await this.get_structure_information();
    await this.click_button('clear');
    await this.wait_for_empty_structure();
    this.replay_tool = null;
    this.replay_atom_indices.clear();
    this.replay_bond_indices.clear();
    await this.build_structure(structure, { sort });
  }

  /** Browser equivalent of Squish get_structure_information() for V3000 state. */
  async get_structure_information() {
    return sketcherCoordinates(parseV3000Structure(await exportText(this.page, 'MDL_MOLV3000')));
  }

  /** Wait only until Qt/WASM has rendered the preceding visible gesture. */
  async wait_for_rendered_object(type, index) {
    for (let attempt = 0; attempt < 100; attempt += 1) {
      const functionName = type === 'atom' ? '_sketcher_get_atom_rect' : '_sketcher_get_bond_rect';
      const rect = await this.page.evaluate(
        ({ name, itemIndex }) => JSON.parse(Module[name](itemIndex)),
        { name: functionName, itemIndex: index },
      );
      if (rect?.width !== undefined) return;
      await this.page.waitForTimeout(10);
    }
    throw new Error(`Timed out waiting for rendered ${type}: ${index}`);
  }

  async wait_for_empty_structure() {
    for (let attempt = 0; attempt < 100; attempt += 1) {
      if (await isEmpty(this.page)) return;
      await this.page.waitForTimeout(10);
    }
    throw new Error('Timed out waiting for Clear to empty the structure');
  }

  /** Match Squish click_tool(): do not re-click an already active sticky tool. */
  async select_replay_tool(tool, widgetName) {
    if (this.replay_tool === tool) return;
    await clickWidget(this.page, widgetName);
    this.replay_tool = tool;
  }

  /** Draw an atom at the coordinate produced by the source Squish wrapper. */
  async add_to_sketcher(x, y, element = 'C', charge = 0, atomCount = null) {
    const tool = TOOL_NAMES[element];
    if (!tool) throw new Error(`Visible structure replay does not yet support element: ${element}`);
    const center = await drawingAreaCenter(this.page);
    const point = { x: center.x + x, y: center.y - y };
    await this.select_replay_tool(element, tool);
    await mouseClick(this.page, point.x, point.y);
    if (atomCount !== null) await this.wait_for_rendered_object('atom', atomCount);
    const chargeTool = charge > 0 ? 'increase_charge_btn' : 'decrease_charge_btn';
    for (let click = 0; click < Math.abs(charge); click += 1) {
      await clickWidget(this.page, chargeTool);
      await mouseClick(this.page, point.x, point.y);
    }
  }

  /** Draw a source V3000 bond with the corresponding visible bond tool. */
  async add_bond(bond, atomsByIndex, bondCount = null) {
    if (bond.order !== 1 && bond.order !== 2) {
      throw new Error(`Visible structure replay does not yet support bond order: ${bond.order}`);
    }
    const first = atomsByIndex.get(bond.atom1);
    const second = atomsByIndex.get(bond.atom2);
    if (!first || !second) throw new Error(`Bond ${bond.index} refers to a missing atom`);
    // bond_order_btn's active/default tool is Double Bond. Opening its popup
    // is only necessary for non-default orders such as Triple Bond.
    await this.select_replay_tool(`bond-${bond.order}`, BOND_TOOL_NAMES[bond.order]);
    const center = await drawingAreaCenter(this.page);
    let start = first;
    let end = second;
    // Exact Molviewer workaround from Squish add_bond (SKETCH-2333).
    if (second.element === 'O') {
      start = { ...second, x: second.x - 20 };
      end = first;
    } else if (first.element === 'O') {
      end = { ...second, x: second.x - 20 };
    }
    await mouseDrag(
      this.page,
      { x: center.x + start.x, y: center.y - start.y },
      { x: center.x + end.x, y: center.y - end.y },
    );
    if (bondCount !== null) await this.wait_for_rendered_object('bond', bondCount);
  }

  /** Browser equivalent of Squish build_structure(..., sort=True/False). */
  async build_structure(structure, { sort = false } = {}) {
    const atoms = [...structure.atoms];
    const bonds = [...structure.bonds];
    if (sort) {
      atoms.sort((first, second) => first.element.localeCompare(second.element));
      bonds.sort((first, second) => first.order - second.order);
    }
    for (let position = 0; position < atoms.length; position += 1) {
      const atom = atoms[position];
      await this.add_to_sketcher(atom.x, atom.y, atom.element, atom.charge, position + 1);
      // Source atom numbers remain stable even when build_structure(sort=True)
      // changes visible creation order.
      this.replay_atom_indices.set(atom.index, position + 1);
    }
    const atomsByIndex = new Map(structure.atoms.map((atom) => [atom.index, atom]));
    for (let position = 0; position < bonds.length; position += 1) {
      await this.add_bond(bonds[position], atomsByIndex, position + 1);
      this.replay_bond_indices.set(bonds[position].index, position + 1);
    }
  }

  /** Equivalent to Squish `click_button(object_name)`. */
  async click_button(object_name) {
    const widgetName = BUTTON_NAMES[object_name] || object_name;
    await clickWidget(this.page, widgetName);
    if (object_name === 'clear') this.current_tool = null;
  }

  /** Equivalent to Squish `click_tool(tool, click_and_hold=False)`. */
  async click_tool(tool, _click_and_hold = false) {
    // Squish tracks sticky tools and avoids re-clicking the current tool.
    // In particular, repeated rect_btn clicks interrupt Shift-add selection.
    if (this.current_tool === tool) return;
    const widgetName = TOOL_NAMES[tool] || tool;
    await clickWidget(this.page, widgetName);
    this.current_tool = tool;
  }

  /** Equivalent to Squish `more_actions_menu(button1, button2=None)`. */
  async more_actions_menu(button1, button2 = null) {
    await clickWidget(this.page, 'more_actions_btn');
    // Qt/WASM creates the popup in a later browser event-loop turn.  Use a
    // Node-side yield here: Playwright's page-side timeout can remain pending
    // while Qt processes a just-triggered model action.
    await new Promise((resolve) => setTimeout(resolve, 100));
    if (button1 === 'modify_all') {
      // Modify All is the first root-menu row.  Popup rows are canvas-only in
      // Qt/WASM, so use the browser-visible popup geometry rather than QMenu.
      // A submenu opens on mouse hover, just as it does for a user.
      await hoverPopupRow(this.page, 0, 18);
      if (button2 === null) {
        return;
      }
      const modifyAllRows = {
        flip_horizontal: 11,
        flip_vertical: 35,
        aromatize: 60,
        kekulize: 83,
        add_explicit_hydrogens: 108,
        remove_explicit_hydrogens: 131,
      };
      const row = modifyAllRows[button2];
      if (row === undefined) throw new Error(`Unsupported Modify All action: ${button2}`);
      await clickPopupRow(this.page, 1, row);
      return;
    }
    if (button1 === 'copy_all_as') {
      await hoverPopupRow(this.page, 0, 202);
      if (button2 === null) {
        return;
      }
      const copyAllAsRows = {
        sdf: 11,
        smi: 35,
        cxsmi: 60,
        inchi: 83,
        inchikey: 108,
        pdb: 131,
      };
      const row = copyAllAsRows[button2];
      if (row === undefined) throw new Error(`Unsupported Copy All As format: ${button2}`);
      await clickPopupRow(this.page, 1, row);
      return;
    }
    if (button1 === 'fit_to_screen') {
      await clickPopupRow(this.page, 0, 249);
      return;
    }
    const rootRows = {
      undo: 35,
      redo: 60,
      select_all: 85,
      clear_selection: 108,
      invert_selection: 131,
      cut: 155,
      copy_all: 179,
      paste: 225,
    };
    if (rootRows[button1] !== undefined) {
      await clickPopupRow(this.page, 0, rootRows[button1]);
      return;
    }
    await clickMenuAction(this.page, MORE_ACTION_NAMES[button1] || button1);
  }

  /** Equivalent to Squish `getClipboardText()` after Copy/Cut actions. */
  async clipboard_text() {
    return clipboardText(this.page);
  }

  /** Equivalent to Squish `type_text("sketcher_area", "<Ctrl+X>")`. */
  async type_text(object_name, text) {
    if (object_name !== 'sketcher_area') {
      throw new Error(`Standalone type_text does not yet support: ${object_name}`);
    }
    await focusCanvas(this.page);
    const shortcut = String(text)
      .replace(/^<Ctrl\+/, 'Control+')
      .replace(/^<Ctrl\+Shift\+/, 'Control+Shift+')
      .replace(/>$/, '');
    await this.page.keyboard.press(shortcut);
  }

  /** Equivalent to `mouse_drag` using source coordinates relative to canvas center. */
  async mouse_drag(x, y, dx, dy, modifier = null, click = 'left') {
    const center = await drawingAreaCenter(this.page);
    const modifiers = modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [];
    await mouseDrag(
      this.page,
      { x: center.x + x, y: center.y + y },
      { x: center.x + x + dx, y: center.y + y + dy },
      { button: click, modifiers },
    );
  }

  /**
   * Return a rendered atom/bond rectangle in the canvas coordinate system.
   * The isolated C++ bridge maps the live QGraphicsItem through QGraphicsView;
   * this avoids the fragile SDF-coordinate-to-pixel calculation in Squish.
   */
  async rendered_object_rect(type, index) {
    const functionName = type === 'atom' ? '_sketcher_get_atom_rect' : '_sketcher_get_bond_rect';
    const rect = await this.page.evaluate(
      ({ name, itemIndex }) => JSON.parse(Module[name](itemIndex)),
      { name: functionName, itemIndex: index },
    );
    if (!rect || rect.width === undefined) {
      throw new Error(`Sketcher ${type} not found: ${index}`);
    }
    return rect;
  }

  /** Equivalent to Squish `click_atom(atom, select=False, modifier=None)`. */
  async click_atom(atom, select = false, modifier = null) {
    if (select) await this.click_tool('rect_btn');
    const renderedIndex = this.replay_atom_indices.get(atom) || atom;
    const rect = await this.rendered_object_rect('atom', renderedIndex);
    await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      modifiers: modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [],
    });
  }

  /** Equivalent to Squish `click_bond(bond, select=False, modifier=None)`. */
  async click_bond(bond, select = false, modifier = null) {
    if (select) await this.click_tool('rect_btn');
    const renderedIndex = this.replay_bond_indices.get(bond) || bond;
    const rect = await this.rendered_object_rect('bond', renderedIndex);
    await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      modifiers: modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [],
    });
  }
}
