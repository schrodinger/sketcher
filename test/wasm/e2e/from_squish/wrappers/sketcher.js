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
  hoverMenuAction,
  isEmpty,
  loadStructureForTest,
  mouseClick,
  mouseDrag,
  openSketcher,
  popupCanvasGeometry,
  sendWidgetMousePress,
  setClipboardText,
  setFilePickerResult,
  widgetRect,
  widgetState,
} from './sketcher_wasm.js';

const BUTTON_NAMES = {
  clear: 'clear_btn',
  clear_selection: 'clear_selection_btn',
  explicit_h: 'explicit_h_btn',
  minus_charge: 'decrease_charge_btn',
  plus_charge: 'increase_charge_btn',
  cleanup: 'cleanup_btn',
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
  erase: 'erase_btn',
  down: 'stereo_bond2_btn',
  single: 'single_bond_btn',
  up: 'stereo_bond1_btn',
  last_picked_element: 'last_picked_element_btn',
  atom_chain: 'atom_chain_btn',
  plus_charge: 'increase_charge_btn',
  minus_charge: 'decrease_charge_btn',
  attachment_point: 'attachment_point_btn',
  cyclohexane: 'cyclohexane_btn',
  benzene: 'benzene_btn',
  cycloheptane: 'cycloheptane_btn',
  cyclopentane: 'cyclopentane_btn',
  cyclopentadiene: 'cyclopentadiene_btn',
  cyclooctane: 'cyclooctane_btn',
  cyclobutane: 'cyclobutane_btn',
  cyclopropane: 'cyclopropane_btn',
};

const BOND_TOOL_NAMES = {
  1: 'single_bond_btn',
  2: 'bond_order_btn',
};

const WILDCARD_TOOL_NAMES = Object.fromEntries(
  ['A', 'AH', 'Q', 'QH', 'M', 'MH', 'X', 'XH'].map((tool) => [
    tool,
    ['atom_query_btn', `${tool}_btn`],
  ]),
);

const POPUP_TOOL_NAMES = {
  single_either: ['stereo_bond2_btn', 'single_either_btn'],
  double_either: ['stereo_bond2_btn', 'double_either_btn'],
  coordinate: ['bond_order_btn', 'coordinate_btn'],
  double: ['bond_order_btn', 'double_btn'],
  triple: ['bond_order_btn', 'triple_btn'],
  zero: ['bond_order_btn', 'zero_btn'],
  aromatic: ['bond_query_btn', 'aromatic_btn'],
  any: ['bond_query_btn', 'any_btn'],
  single_double: ['bond_query_btn', 'single_double_btn'],
  single_aromatic: ['bond_query_btn', 'single_aromatic_btn'],
  double_aromatic: ['bond_query_btn', 'double_aromatic_btn'],
  Pd: ['periodic_table_btn', 'pd_btn'],
};

// The periodic-table widget follows the same object-name convention for every
// element. Keep the complete Squish tool list here so `tst_tools` and hidden
// shortcuts use the source `send_event()` -> real child-button interaction.
const PERIODIC_TABLE_TOOLS = [
  'He',
  'Li',
  'Be',
  'B',
  'Ne',
  'Na',
  'Mg',
  'Al',
  'Si',
  'Ar',
  'K',
  'Ca',
  'Sc',
  'Ti',
  'V',
  'Cr',
  'Mn',
  'Fe',
  'Co',
  'Ni',
  'Cu',
  'Zn',
  'Ga',
  'Ge',
  'As',
  'Se',
  'Br',
  'Kr',
  'Rb',
  'Sr',
  'Y',
  'Zr',
  'Nb',
  'Mo',
  'Tc',
  'Ru',
  'Rh',
  'Pd',
  'Ag',
  'Cd',
  'In',
  'Sn',
  'Sb',
  'Te',
  'I',
  'Xe',
  'Cs',
  'Ba',
  'La',
  'Ce',
  'Pr',
  'Nd',
  'Pm',
  'Sm',
  'Eu',
  'Gd',
  'Tb',
  'Dy',
  'Ho',
  'Er',
  'Tm',
  'Yb',
  'Lu',
  'Hf',
  'Ta',
  'W',
  'Re',
  'Os',
  'Ir',
  'Pt',
  'Au',
  'Hg',
  'Tl',
  'Pb',
  'Bi',
  'Po',
  'At',
  'Rn',
  'Fr',
  'Ra',
  'Ac',
  'Th',
  'Pa',
  'U',
  'Np',
  'Pu',
  'Am',
  'Cm',
  'Bk',
  'Cf',
  'Es',
  'Fm',
  'Md',
  'No',
  'Lr',
  'Rf',
  'Db',
  'Sg',
  'Bh',
  'Hs',
  'Mt',
  'Ds',
  'Rg',
  'Cn',
  'Nh',
  'Fl',
  'Mc',
  'Lv',
  'Ts',
  'Og',
];
for (const element of PERIODIC_TABLE_TOOLS) {
  POPUP_TOOL_NAMES[element] = ['periodic_table_btn', `${element.toLowerCase()}_btn`];
}

const RGROUP_TOOLS = ['r+', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9'];
for (const rgroup of RGROUP_TOOLS) {
  const number = rgroup === 'r+' ? 'new' : rgroup.slice(1);
  POPUP_TOOL_NAMES[rgroup] = [
    'rgroup_btn',
    `${number === 'new' ? 'new_rgroup' : `rgroup${number}`}_btn`,
  ];
}

const REACTION_TOOLS = ['rxn_arrow', 'rxn_plus', 'add_mapping', 'remove_mapping'];
for (const reactionTool of REACTION_TOOLS) {
  POPUP_TOOL_NAMES[reactionTool] = ['reaction_btn', `${reactionTool}_btn`];
}

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

// Squish passes stable snake-case identifiers; context-menu actions expose
// their human-visible labels.  Keep this translation in the wrapper so test
// bodies retain the source call format.
const CONTEXT_MENU_NAMES = {
  copy_as: 'Copy As',
  modify_atoms: 'Modify Atoms',
  modify_bonds: 'Modify Bonds',
  replace_atoms_with: 'Replace Atoms with',
  set_element: 'Set Element',
  wildcard: 'Wildcard',
  other_type: 'Other Type',
  single_up_down: 'Single Up/Down',
  double_cis_trans: 'Double Cis/Trans',
  zero_order: 'Zero Order',
  'Not In a Ring': 'Not In a Ring',
  '+_charge': '+ Charge',
  '–_charge': '– Charge',
  add_explicit_hydrogens: 'Add Explicit Hydrogens',
  remove_explicit_hydrogens: 'Remove Explicit Hydrogens',
  add_unpaired_electron: 'Add Unpaired Electron',
  remove_unpaired_electron: 'Remove Unpaired Electron',
  A: 'A (Any heavy atom)',
  AH: 'AH (Any or H)',
  Q: 'Q (Heteroatom)',
  QH: 'QH (Hetero or H)',
  M: 'M (Metal)',
  MH: 'MH (Metal or H)',
  X: 'X (Halogen)',
  XH: 'XH (Halogen or H)',
};

const IMPORT_MENU_NAMES = {
  import_from_file: 'Import from File...',
  paste_in_text: 'Paste in Text...',
  replace_current_content: 'Replace Current Content',
};

/** Return the largest visible Qt popup canvas: the active modal dialog. */
async function dialogCanvasGeometry(page) {
  for (let attempt = 0; attempt < 50; attempt += 1) {
    const dialogs = await page.locator('[id^="qt-window-"] canvas').evaluateAll((canvases) =>
      canvases
        .map((canvas) => {
          const rect = canvas.getBoundingClientRect();
          return { height: rect.height, width: rect.width, x: rect.x, y: rect.y };
        })
        .filter((rect) => rect.width < window.innerWidth || rect.height < window.innerHeight)
        .sort((first, second) => second.width * second.height - first.width * first.height),
    );
    if (dialogs[0]?.width && dialogs[0]?.height) return dialogs[0];
    await page.waitForTimeout(25);
  }
  throw new Error('Paste in Text dialog did not expose browser canvas geometry');
}

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
    this.replay_atom_coordinates = new Map();
    this.replay_bond_coordinates = new Map();
    // Affine conversion from Squish's centred drawing coordinates to the
    // current browser canvas.  It is learned from the visible imported atoms,
    // so tests can retain their source gestures without a clear/rebuild step.
    this.replay_canvas_transform = null;
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
    this.replay_atom_coordinates.clear();
    this.replay_bond_coordinates.clear();
    this.replay_canvas_transform = null;
    await this.build_structure(structure, { sort });
  }

  /** Browser equivalent of Squish get_structure_information() for V3000 state. */
  async get_structure_information() {
    // MolViewer's Copy All returns V3000; standalone Sketcher returns SMILES.
    // Use the visible Copy All As -> SDF menu path to obtain the same V3000
    // data needed for the source's clear-and-visible-rebuild workflow.
    // Clear any prior browser clipboard payload first. Otherwise a previous
    // V3000 export can satisfy the read loop before Qt/WASM completes this
    // user-visible Copy As action.
    await setClipboardText(this.page, '__playwright_waiting_for_sdf__');
    await this.more_actions_menu('copy_all_as', 'sdf');
    let structureText = '';
    for (let attempt = 0; attempt < 50; attempt += 1) {
      structureText = await this.clipboard_text();
      if (structureText.includes('M  V30 BEGIN ATOM')) break;
      await this.page.waitForTimeout(25);
    }
    if (!structureText.includes('M  V30 BEGIN ATOM')) {
      throw new Error(
        `More Actions -> Copy All As SDF returned unexpected clipboard data: ${structureText}`,
      );
    }
    return sketcherCoordinates(parseV3000Structure(structureText));
  }

  /** Map V3000 source atom numbers to the live imported RDKit atom indexes. */
  async map_imported_atom_indexes() {
    const sdf = await this.copy_all_as_text('sdf');
    const source = parseV3000Structure(sdf);
    const sourceCoordinates = sketcherCoordinates(source);
    const rendered = await this.page.evaluate(() =>
      JSON.parse(Module._sketcher_get_rendered_atom_geometry()),
    );
    const renderedBonds = await this.page.evaluate(() =>
      JSON.parse(Module._sketcher_get_rendered_bond_geometry()),
    );
    this.replay_atom_indices.clear();
    this.replay_atom_coordinates.clear();
    const transformPoints = [];
    for (const atom of source.atoms) {
      const match = rendered.find(
        (candidate) =>
          Math.abs(candidate.x - atom.x) < 1e-5 && Math.abs(candidate.y - atom.y) < 1e-5,
      );
      if (!match)
        throw new Error(`Unable to map imported V3000 atom ${atom.index} to live geometry`);
      this.replay_atom_indices.set(atom.index, match.index);
      // Retain the exact coordinate system returned by the Squish wrapper's
      // get_structure_information(), then map it to the live browser canvas.
      // This avoids changing the source gesture merely because WASM lays out
      // the imported molecule at a different pixel location.
      const squishAtom = sourceCoordinates.atoms.find(
        (candidate) => candidate.index === atom.index,
      );
      this.replay_atom_coordinates.set(atom.index, { x: squishAtom.x, y: squishAtom.y });
      const rect = await this.rendered_object_rect('atom', match.index);
      transformPoints.push({
        sourceX: squishAtom.x,
        sourceY: squishAtom.y,
        canvasX: rect.x + rect.width / 2,
        canvasY: rect.y + rect.height / 2,
      });
    }
    this.replay_canvas_transform = this.canvas_transform_from_points(transformPoints);
    this.replay_bond_indices.clear();
    for (const bond of source.bonds) {
      const atom1 = this.replay_atom_indices.get(bond.atom1);
      const atom2 = this.replay_atom_indices.get(bond.atom2);
      const match = renderedBonds.find(
        (candidate) =>
          (candidate.atom1 === atom1 && candidate.atom2 === atom2) ||
          (candidate.atom1 === atom2 && candidate.atom2 === atom1),
      );
      if (!match) {
        throw new Error(`Unable to map imported V3000 bond ${bond.index} to live geometry`);
      }
      this.replay_bond_indices.set(bond.index, match.index);
    }
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
      this.replay_atom_coordinates.set(atom.index, { x: atom.x, y: atom.y });
    }
    const atomsByIndex = new Map(structure.atoms.map((atom) => [atom.index, atom]));
    for (let position = 0; position < bonds.length; position += 1) {
      await this.add_bond(bonds[position], atomsByIndex, position + 1);
      this.replay_bond_indices.set(bonds[position].index, position + 1);
      const first = atomsByIndex.get(bonds[position].atom1);
      const second = atomsByIndex.get(bonds[position].atom2);
      this.replay_bond_coordinates.set(bonds[position].index, {
        x: (first.x + second.x) / 2,
        y: (first.y + second.y) / 2,
      });
    }
  }

  /** Equivalent to Squish `click_button(object_name)`. */
  async click_button(object_name) {
    if (object_name === 'about_sketcher_ok') {
      // About2DSketcher is a Qt/WASM top-level canvas. Its generated
      // QDialogButtonBox is not consistently discoverable by objectName after
      // the Help menu's queued action, but its rendered, centered OK button is
      // stably positioned at the bottom of that live dialog canvas.
      const dialog = await popupCanvasGeometry(this.page, 0);
      await mouseClick(this.page, dialog.x + dialog.width / 2, dialog.y + dialog.height - 16);
      return;
    }
    const widgetName = BUTTON_NAMES[object_name] || object_name;
    await clickWidget(this.page, widgetName);
    if (object_name === 'clear') this.current_tool = null;
  }

  /** Equivalent to Squish `click_tool(tool, click_and_hold=False)`. */
  async click_tool(tool, click_and_hold = false) {
    // Squish tracks sticky tools and avoids re-clicking the current tool.
    // In particular, repeated rect_btn clicks interrupt Shift-add selection.
    if (this.current_tool === tool) return;
    if (WILDCARD_TOOL_NAMES[tool]) {
      const [parent, child] = WILDCARD_TOOL_NAMES[tool];
      if (!click_and_hold) {
        await clickWidget(this.page, parent);
        await clickWidget(this.page, parent);
        await this.page.waitForTimeout(100);
        try {
          await clickWidget(this.page, child);
        } catch (error) {
          // Qt/WASM can discard the just-opened popup geometry for one frame.
          // Reopen it with the same ordinary human click rather than invoking
          // the child action directly.
          await clickWidget(this.page, parent);
          await this.page.waitForTimeout(50);
          await clickWidget(this.page, child);
        }
        this.current_tool = tool;
        return;
      }
      const rect = await widgetRect(this.page, parent);
      await this.page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2, { steps: 4 });
      await this.page.mouse.down();
      await this.page.waitForTimeout(600);
      try {
        await clickWidget(this.page, child);
      } finally {
        await this.page.mouse.up();
      }
      this.current_tool = tool;
      return;
    }
    if (POPUP_TOOL_NAMES[tool]) {
      const [parent, child] = POPUP_TOOL_NAMES[tool];
      if (PERIODIC_TABLE_TOOLS.includes(tool)) {
        // Exact source path: Sketcher.click_tool() opens the periodic-table
        // ToolButtonWithPopup with Squish send_event(), then clicks the
        // element child.  A normal browser click selects the parent tool
        // instead, so retain this isolated test-only Qt press shim.
        await sendWidgetMousePress(this.page, parent);
        await this.page.waitForTimeout(25);
        await clickWidget(this.page, child);
        this.current_tool = tool;
        return;
      }
      if (!click_and_hold) {
        // Source click_tool(..., click_and_hold=False) uses the ordinary
        // two-click ToolButtonWithPopup path: choose the parent, reopen its
        // popup, then choose the child.  Preserve those browser clicks rather
        // than silently substituting a direct Qt action.
        await clickWidget(this.page, parent);
        await clickWidget(this.page, parent);
        await this.page.waitForTimeout(100);
        try {
          await clickWidget(this.page, child);
        } catch (error) {
          // Reopen a transiently unavailable Qt/WASM popup with the same
          // browser click sequence used by the source wrapper.
          await clickWidget(this.page, parent);
          await this.page.waitForTimeout(50);
          await clickWidget(this.page, child);
        }
        this.current_tool = tool;
        return;
      }
      const rect = await widgetRect(this.page, parent);
      await this.page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2, { steps: 4 });
      await this.page.mouse.down();
      // ToolButtonWithPopup uses a 250 ms hold timer.  Keep a margin so the
      // browser event loop can deliver the same popup-opening behavior that
      // Squish's press-only sendEvent relies on.
      await this.page.waitForTimeout(600);
      try {
        await clickWidget(this.page, child);
      } finally {
        await this.page.mouse.up();
      }
      this.current_tool = tool;
      return;
    }
    if (tool === 'lasso_btn' || tool === 'ellipse_btn') {
      // A checked ToolButtonWithPopup opens its visible selector popup on a
      // second ordinary click.  Choose the real popup child with the mouse.
      if (
        this.current_tool !== 'rect_btn' &&
        this.current_tool !== 'lasso_btn' &&
        this.current_tool !== 'ellipse_btn'
      ) {
        await clickWidget(this.page, 'select_tool_btn');
      }
      await clickWidget(this.page, 'select_tool_btn');
      await this.page.waitForTimeout(100);
      await clickWidget(this.page, tool);
      this.current_tool = tool;
      return;
    }
    const widgetName = TOOL_NAMES[tool] || tool;
    await clickWidget(this.page, widgetName);
    this.current_tool = tool;
  }

  /** Return the visible Qt button Squish observes for a sticky tool. */
  tool_widget_name(tool) {
    if (PERIODIC_TABLE_TOOLS.includes(tool)) return `${tool.toLowerCase()}_btn`;
    if (WILDCARD_TOOL_NAMES[tool]) return WILDCARD_TOOL_NAMES[tool][0];
    if (POPUP_TOOL_NAMES[tool]) return POPUP_TOOL_NAMES[tool][0];
    return TOOL_NAMES[tool] || tool;
  }

  /** Read the visible Qt state Squish verifies for a mapped tool/button. */
  async widget_state(object_name) {
    return widgetState(this.page, object_name);
  }

  /**
   * Source-compatible reset between independent Squish test sections.
   *
   * The standalone WASM app has the same default view options as the source
   * MolViewer session, so this preserves the source's tool reset and clear
   * sequence while leaving the already-default configuration untouched.
   */
  async reset_state() {
    await this.click_button('clear_selection');
    await this.click_tool('rect_btn');
    await this.click_tool('A');
    await this.click_tool('Si');
    await this.click_tool('down');
    await this.click_tool('double');
    await this.click_tool('aromatic');
    await this.click_tool('r+');
    await this.click_tool('rxn_arrow');
    await this.click_tool('C');
    await this.click_button('clear');
    this.replay_atom_indices.clear();
    this.replay_bond_indices.clear();
    this.replay_atom_coordinates.clear();
    this.replay_bond_coordinates.clear();
    this.replay_canvas_transform = null;
  }

  /** Move the real browser mouse to a source-coordinate canvas point. */
  async mouse_move(x = 0, y = 0) {
    const point = await this.source_canvas_point(x, y);
    await this.page.mouse.move(point.x, point.y, { steps: 4 });
  }

  /** Source-compatible Configure View menu entry point. */
  async configure_view_menu(button, options = {}) {
    const labels = {
      valence_errors: 'Valence Errors',
      heteroatom_colors: 'Heteroatom Colors',
      stereo_labels: 'Stereo Labels',
      preferences: 'Preferences',
    };
    await clickWidget(this.page, 'configure_view_btn');
    await this.page.waitForTimeout(100);
    await clickMenuAction(this.page, labels[button] || button);
    if (button !== 'preferences') return;
    const {
      atom_font_size = 18,
      bond_line_width = 2,
      label_carbons = false,
      terminal_only = true,
      all_carbons = false,
      show_stereo_annotations = true,
      with_abs = false,
      include_undefined_centers = true,
      color_heteroatoms = true,
      color_mode = 'default',
      reset = false,
      close = true,
    } = options;
    if (atom_font_size !== 18)
      await setWidgetText(this.page, 'm_atom_font_size_sb', atom_font_size);
    if (bond_line_width !== 2)
      await setWidgetText(this.page, 'm_bond_line_width_sb', bond_line_width);
    if (label_carbons) {
      await clickWidget(this.page, 'm_label_carbons_cb');
      if (all_carbons === terminal_only)
        throw new Error('Choose exactly one carbon-label radio option');
      if (all_carbons) await clickWidget(this.page, 'm_label_all_C_rb');
    }
    if (!show_stereo_annotations) await clickWidget(this.page, 'm_show_stereo_cb');
    if (with_abs) await clickWidget(this.page, 'm_abs_cb');
    if (!include_undefined_centers) await clickWidget(this.page, 'm_undefined_centers_labels_cb');
    if (!color_heteroatoms) await clickWidget(this.page, 'm_color_heteroatoms_cb');
    if (color_heteroatoms && color_mode.toLowerCase() !== 'default') {
      await clickWidget(this.page, 'm_color_mode_combo');
      await this.page.keyboard.type(String(color_mode));
      await this.page.keyboard.press('Enter');
    }
    if (reset) await clickWidget(this.page, 'm_reset_to_default_btn');
    if (close) await clickWidget(this.page, 'm_close_btn');
  }

  /** Equivalent to Squish `help_menu(button)`. */
  async help_menu(button) {
    const rows = {
      getting_started: 35,
      about_sketcher: 58,
    };
    if (rows[button] === undefined) throw new Error(`Unsupported Help menu action: ${button}`);
    await clickWidget(this.page, 'sketcher_help_btn');
    await this.page.waitForTimeout(100);
    // Qt/WASM blocks when its QMenu action geometry is queried during this
    // popup's lifecycle. Use its browser-visible popup canvas, as a human
    // pointer does, rather than invoking the menu action directly.
    await clickPopupRow(this.page, 0, rows[button]);
  }

  /** First standalone slice of Squish `edit_atom_properties()`. */
  async edit_atom_properties({
    set_as = 'atom',
    query_type,
    dropdown_type,
    element,
    isotope,
    charge,
    unpaired_electrons,
    enhanced_stereo_label,
    element_list,
    click_ok = true,
    click_cancel = false,
    click_reset = false,
  } = {}) {
    if (set_as === 'query') await clickWidget(this.page, 'set_as_query_rb');
    else await clickWidget(this.page, 'set_as_atom_rb');
    for (const [widget, value] of [
      ['query_type_combo', query_type],
      ['query_type_combo', dropdown_type],
    ]) {
      if (value === undefined) continue;
      await clickWidget(this.page, widget);
      await this.page.keyboard.type(String(value));
      await this.page.keyboard.press('Enter');
    }
    if (element !== undefined) await setWidgetText(this.page, 'element_le', element);
    if (isotope !== undefined) await setWidgetText(this.page, 'isotope_le', isotope);
    if (charge !== undefined) await setWidgetText(this.page, 'charge_sb', charge);
    if (unpaired_electrons !== undefined)
      await setWidgetText(this.page, 'unpaired_elec_sb', unpaired_electrons);
    if (enhanced_stereo_label !== undefined)
      await setWidgetText(this.page, 'enhanced_stereo_sb', enhanced_stereo_label);
    if (element_list !== undefined) await setWidgetText(this.page, 'element_list_le', element_list);
    if (click_reset) await clickWidget(this.page, 'Reset');
    if (click_cancel) await clickWidget(this.page, 'Cancel');
    else if (click_ok) await clickWidget(this.page, 'OK');
  }

  /**
   * Browser equivalent of Squish `import_menu()` for the standalone Sketcher.
   * Menu activation, text entry, and confirmation are all real user input.
   */
  async import_menu(button, text = null, close_panel = true) {
    const menuItem = IMPORT_MENU_NAMES[button] || button;
    if (button === 'import_from_file' && text !== null) {
      await setFilePickerResult(this.page, String(text));
    }
    const importRows = {
      import_from_file: 12,
      paste_in_text: 35,
      replace_current_content: 58,
    };
    if (importRows[button] !== undefined) {
      let lastError;
      for (let attempt = 0; attempt < 3; attempt += 1) {
        await clickWidget(this.page, 'import_btn');
        await new Promise((resolve) => setTimeout(resolve, 100));
        try {
          await clickPopupRow(this.page, 0, importRows[button]);
          lastError = null;
          break;
        } catch (error) {
          lastError = error;
          await this.page.waitForTimeout(50);
        }
      }
      if (lastError) throw lastError;
    } else {
      await clickWidget(this.page, 'import_btn');
      await new Promise((resolve) => setTimeout(resolve, 100));
      await clickMenuAction(this.page, menuItem);
    }

    if (button === 'replace_current_content') {
      // The checkable menu action updates the dialog's import mode on the
      // next Qt event-loop turn.  Let that user-visible menu action settle
      // before opening Import again.
      await this.page.waitForTimeout(100);
    }

    if (button === 'import_from_file' && text !== null) {
      await this.page.waitForTimeout(100);
      return;
    }

    if (button !== 'paste_in_text' || text === null) {
      return;
    }

    // Match the Squish workflow: seed the system clipboard, then paste through
    // the dialog's focused text editor with an actual Ctrl+V key gesture.
    await setClipboardText(this.page, String(text));
    await clickWidget(this.page, 'structure_text_edit');
    await this.page.keyboard.press('Control+V');
    if (close_panel) {
      await this.click_dialog_ok();
    }
  }

  /** Return the displayed Paste in Text dialog status, not model state. */
  async paste_in_text_status() {
    return (await widgetState(this.page, 'status_lbl')).text;
  }

  /** Return the visible autodetected format label from Paste in Text. */
  async paste_in_text_format() {
    return (await widgetState(this.page, 'format_combo')).text;
  }

  /** Click the standard OK half of a visible Qt QDialogButtonBox. */
  async click_dialog_ok() {
    const dialog = await dialogCanvasGeometry(this.page);
    // Visible Qt/WASM dialog layout: OK precedes Cancel in the bottom-right
    // button row.  Use the dialog canvas's browser geometry, not Qt internals.
    await mouseClick(this.page, dialog.x + dialog.width * 0.71, dialog.y + dialog.height * 0.96);
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
        // The standalone menu includes Maestro between SDF and SMILES.
        // Keep these coordinates aligned with get_standard_export_formats().
        smi: 60,
        cxsmi: 83,
        inchi: 155,
        inchikey: 179,
        pdb: 202,
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

  /**
   * Activate the visible Copy All As action and wait for its clipboard result.
   *
   * Qt/WASM updates the browser clipboard asynchronously.  A sentinel avoids
   * accidentally treating a preceding Ctrl+C payload as this action's result.
   */
  async copy_all_as_text(format) {
    const sentinel = `__playwright_waiting_for_${format}__`;
    await setClipboardText(this.page, sentinel);
    await this.more_actions_menu('copy_all_as', format);
    return this.wait_for_clipboard_change(sentinel, `Copy All As ${format}`);
  }

  /** Equivalent to context-menu Copy As followed by Squish getClipboardText(). */
  async selection_context_copy_as_text(target, format) {
    const sentinel = `__playwright_waiting_for_selection_${format}__`;
    await setClipboardText(this.page, sentinel);
    await this.selection_context_menu(target, 'copy_as', format);
    return this.wait_for_clipboard_change(sentinel, `Selection Copy As ${format}`);
  }

  /** Equivalent to context-menu Copy followed by Squish getClipboardText(). */
  async selection_context_copy_text(target) {
    const sentinel = '__playwright_waiting_for_selection_copy__';
    await setClipboardText(this.page, sentinel);
    await this.selection_context_menu(target, 'copy');
    return this.wait_for_clipboard_change(sentinel, 'Selection Copy');
  }

  /** Wait for Qt/WASM to replace a browser-clipboard sentinel. */
  async wait_for_clipboard_change(sentinel, action = 'Sketcher action') {
    for (let attempt = 0; attempt < 50; attempt += 1) {
      const text = await this.clipboard_text();
      if (text !== sentinel) return text;
      await this.page.waitForTimeout(25);
    }
    throw new Error(`${action} did not update the browser clipboard`);
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
    const sourceKey = String(text);
    const shortcut = sourceKey
      .replace(/^<Ctrl\+/, 'Control+')
      .replace(/^<Ctrl\+Shift\+/, 'Control+Shift+')
      .replace(/>$/, '');
    const key =
      {
        '<Backspace>': 'Backspace',
        ' ': 'Space',
      }[sourceKey] || shortcut;
    await focusCanvas(this.page);
    await this.page.keyboard.press(key);
  }

  /** Equivalent to `mouse_drag` using source coordinates relative to canvas center. */
  async mouse_drag(x, y, dx, dy, modifier = null, click = 'left') {
    const modifiers = modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [];
    await mouseDrag(
      this.page,
      await this.source_canvas_point(x, y),
      await this.source_canvas_point(x + dx, y + dy),
      { button: click, modifiers },
    );
  }

  /** Learn an axis-aligned source-coordinate to browser-coordinate mapping. */
  canvas_transform_from_points(points) {
    const fit = (sourceKey, canvasKey) => {
      const sourceMean = points.reduce((sum, point) => sum + point[sourceKey], 0) / points.length;
      const canvasMean = points.reduce((sum, point) => sum + point[canvasKey], 0) / points.length;
      const variance = points.reduce((sum, point) => sum + (point[sourceKey] - sourceMean) ** 2, 0);
      const covariance = points.reduce(
        (sum, point) => sum + (point[sourceKey] - sourceMean) * (point[canvasKey] - canvasMean),
        0,
      );
      const scale = variance ? covariance / variance : 1;
      return { scale, offset: canvasMean - scale * sourceMean };
    };
    return { x: fit('sourceX', 'canvasX'), y: fit('sourceY', 'canvasY') };
  }

  /** Convert a source Squish drawing-area coordinate to its browser pixel. */
  async source_canvas_point(x, y) {
    if (this.replay_canvas_transform) {
      const view = await widgetRect(this.page, 'view');
      const point = {
        x: this.replay_canvas_transform.x.offset + this.replay_canvas_transform.x.scale * x,
        y: this.replay_canvas_transform.y.offset + this.replay_canvas_transform.y.scale * y,
      };
      return {
        x: Math.max(view.x + 2, Math.min(view.x + view.width - 2, point.x)),
        y: Math.max(view.y + 2, Math.min(view.y + view.height - 2, point.y)),
      };
    }
    const center = await drawingAreaCenter(this.page);
    return { x: center.x + x, y: center.y - y };
  }

  /**
   * Equivalent to Squish `click_sketcher(x, y, select=False, click='left')`.
   *
   * Squish coordinates are relative to the drawing area's center and use an
   * upward-positive Y axis.  The browser canvas can be smaller than the
   * desktop MolViewer canvas, so retain the source intent (a background
   * click) while keeping an out-of-range point inside its visible viewport.
   */
  async click_sketcher(x, y, select = false, click = 'left') {
    if (select) await this.click_tool('rect_btn');
    const rect = await widgetRect(this.page, 'view');
    const point = {
      x: Math.max(rect.x + 2, Math.min(rect.x + rect.width - 2, rect.x + rect.width / 2 + x)),
      y: Math.max(rect.y + 2, Math.min(rect.y + rect.height - 2, rect.y + rect.height / 2 - y)),
    };
    await mouseClick(this.page, point.x, point.y, { button: click });
  }

  /**
   * Activate a top-level selection context-menu action with a real right
   * click followed by a real menu-item click.
   *
   * `target` may be a source-style coordinate pair or `{ type, index }`.
   * The latter resolves the current rendered QGraphicsItem and is used for
   * imported structures, whose desktop Squish coordinates were fragile.
   */
  async selection_context_menu(target, yOrAction, action = null, ...remainingActions) {
    let actionNames;
    let openContextMenu;
    if (typeof target === 'object') {
      actionNames = [yOrAction, action, ...remainingActions].filter((name) => name !== null);
      const rect = await this.rendered_object_rect(
        target.type,
        this.replay_atom_indices.get(target.index) ||
          this.replay_bond_indices.get(target.index) ||
          target.index,
      );
      openContextMenu = () =>
        mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
          button: 'right',
        });
    } else {
      openContextMenu = () => this.click_sketcher(target, yOrAction, false, 'right');
      actionNames = [action, ...remainingActions].filter((name) => name !== null);
    }
    if (!actionNames.length) throw new Error('selection_context_menu requires an action');
    const labelFor = (name) =>
      CONTEXT_MENU_NAMES[name] ||
      COPY_ALL_AS_NAMES[name] ||
      String(name)
        .replace(/_/g, ' ')
        .replace(/\b\w/g, (letter) => letter.toUpperCase());
    let lastError;
    // Qt/WASM occasionally drops a just-opened popup from the bridge for one
    // frame. Reopen the same human context menu rather than falling back to a
    // direct action call; this keeps the operation user-equivalent.
    for (let attempt = 0; attempt < 3; attempt += 1) {
      await openContextMenu();
      try {
        for (const actionName of actionNames.slice(0, -1)) {
          await hoverMenuAction(this.page, labelFor(actionName));
        }
        await clickMenuAction(this.page, labelFor(actionNames.at(-1)));
        return;
      } catch (error) {
        lastError = error;
        await this.page.waitForTimeout(50);
      }
    }
    throw lastError;
  }

  /** Open a selection context submenu and leave its widget content visible. */
  async open_selection_context_submenu(target, ...menuNames) {
    if (!menuNames.length) throw new Error('open_selection_context_submenu requires a menu path');
    const rect = await this.rendered_object_rect(
      target.type,
      this.replay_atom_indices.get(target.index) ||
        this.replay_bond_indices.get(target.index) ||
        target.index,
    );
    const labelFor = (name) =>
      CONTEXT_MENU_NAMES[name] ||
      COPY_ALL_AS_NAMES[name] ||
      String(name)
        .replace(/_/g, ' ')
        .replace(/\b\w/g, (letter) => letter.toUpperCase());
    let lastError;
    for (let attempt = 0; attempt < 3; attempt += 1) {
      await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
        button: 'right',
      });
      try {
        for (const menuName of menuNames) {
          await hoverMenuAction(this.page, labelFor(menuName));
        }
        return;
      } catch (error) {
        lastError = error;
        await this.page.waitForTimeout(50);
      }
    }
    throw lastError;
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

  /** Preserve source object targets before edits renumber live RDKit items. */
  async capture_rendered_targets(type, indexes) {
    const targets = new Map();
    for (const index of indexes) {
      targets.set(index, await this.rendered_object_rect(type, index));
    }
    return targets;
  }

  /** Click a previously captured visible item target with the real mouse. */
  async click_rendered_target(rect) {
    await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2);
  }

  /** Drag from a captured target, clamped to the visible drawing viewport. */
  async drag_from_rendered_target(rect, dx, dy) {
    const view = await widgetRect(this.page, 'view');
    const start = { x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 };
    const end = {
      x: Math.max(view.x + 2, Math.min(view.x + view.width - 2, start.x + dx)),
      y: Math.max(view.y + 2, Math.min(view.y + view.height - 2, start.y + dy)),
    };
    await mouseDrag(this.page, start, end);
  }

  /**
   * Replay a Squish drag that starts relative to a live item and ends at a
   * source-coordinate point in the drawing area.  This is useful for the
   * Move/Rotate circle, whose source gesture intentionally begins off atom.
   */
  async drag_rendered_target_to_source_point(rect, startOffsetX, startOffsetY, targetX, targetY) {
    const start = await this.source_canvas_point(
      this.replay_atom_coordinates.get(18)?.x + startOffsetX,
      this.replay_atom_coordinates.get(18)?.y + startOffsetY,
    );
    await mouseDrag(this.page, start, await this.source_canvas_point(targetX, targetY));
  }

  /** Equivalent to Squish startDrag/dropOn using canvas-local pixels. */
  async drag_canvas_point(startX, startY, endX, endY) {
    const view = await widgetRect(this.page, 'view');
    await mouseDrag(
      this.page,
      { x: view.x + startX, y: view.y + startY },
      { x: view.x + endX, y: view.y + endY },
    );
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
