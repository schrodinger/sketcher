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
  clickWidget,
  drawingAreaCenter,
  focusCanvas,
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
  move_rotate: 'move_rotate_btn',
  rect_btn: 'select_tool_btn',
};

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
  }

  async open() {
    await openSketcher(this.page);
  }

  /** Fixture setup only; never use this method in `tst_import_menu`. */
  async load_structure_for_test(text) {
    await loadStructureForTest(this.page, text);
  }

  /** Equivalent to Squish `click_button(object_name)`. */
  async click_button(object_name) {
    const widgetName = BUTTON_NAMES[object_name] || object_name;
    await clickWidget(this.page, widgetName);
  }

  /** Equivalent to Squish `click_tool(tool, click_and_hold=False)`. */
  async click_tool(tool, _click_and_hold = false) {
    const widgetName = TOOL_NAMES[tool] || tool;
    await clickWidget(this.page, widgetName);
  }

  /** Equivalent to Squish `more_actions_menu(button1, button2=None)`. */
  async more_actions_menu(button1, button2 = null) {
    await clickWidget(this.page, 'more_actions_btn');
    if (button1 === 'modify_all') {
      await clickMenuAction(this.page, 'Modify All');
      if (button2 === null) {
        return;
      }
      await clickMenuAction(this.page, MORE_ACTION_NAMES[button2] || button2);
      return;
    }
    if (button1 === 'copy_all_as') {
      await clickMenuAction(this.page, MORE_ACTION_NAMES[button1]);
      if (button2 === null) {
        return;
      }
      await clickMenuAction(this.page, COPY_ALL_AS_NAMES[button2] || button2);
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
    const rect = await this.rendered_object_rect('atom', atom);
    await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      modifiers: modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [],
    });
  }

  /** Equivalent to Squish `click_bond(bond, select=False, modifier=None)`. */
  async click_bond(bond, select = false, modifier = null) {
    if (select) await this.click_tool('rect_btn');
    const rect = await this.rendered_object_rect('bond', bond);
    await mouseClick(this.page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      modifiers: modifier === 'shift' ? ['Shift'] : modifier === 'control' ? ['Control'] : [],
    });
  }
}
