import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { clickWidget, closeActiveQtPopups, hideMouseMarker, isEmpty, sendWidgetMousePress, widgetState } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
test.setTimeout(180_000);

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  // Qt/WASM may retain a closing context-menu canvas for one render frame.
  await page.waitForTimeout(150);
  await expect(page.locator('#screen')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_select_mode_active_selection', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Move or rotate selected structure — source checkpoint order preserved.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.click_button('select_all');
    await expect((await widgetState(page, 'move_rotate_btn')).enabled).toBe(true);
    await sk.click_button('move_rotate');
    await sk.mouse_drag(0, 0, 0, 50);
    await checkpoint(page, 'simple_move');
    await sk.mouse_drag(-50, -50, 0, 450, null, 'middle');
    await checkpoint(page, 'simple_rotate');

    // Clear Selection — More Actions, Select Options, then a background click.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.click_button('select_all');
    await sk.more_actions_menu('clear_selection');
    await checkpoint(page, 'clear_selection_more_actions_menu');

    await sk.click_button('select_all');
    await sk.click_button('clear_selection');
    await checkpoint(page, 'clear_selection_select_options');

    await sk.click_button('select_all');
    await sk.click_tool('rect_btn');
    await sk.click_sketcher(-500, -500);
    await checkpoint(page, 'clicking_background_clears_selection');

    // Invert Selection — context menu, then More Actions.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    for (const n of [1, 4, 7]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    await sk.selection_context_menu({ type: 'atom', index: 1 }, 'invert_selection');
    await checkpoint(page, 'invert_selection_context');
    await sk.more_actions_menu('invert_selection');
    await checkpoint(page, 'invert_selection_more_actions_menu');

    // Delete Selected Structure — context menu, then the Select Options eraser.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    await sk.click_button('select_all');
    await sk.selection_context_menu({ type: 'atom', index: 1 }, 'delete');
    await expect.poll(() => isEmpty(page)).toBe(true);

    await sk.import_menu('paste_in_text', SOURCE);
    await sk.click_button('select_all');
    await sk.click_tool('erase');
    await expect.poll(() => isEmpty(page)).toBe(true);

    // Edit Selected Structure via Edit Actions: the direct toolbar actions.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    for (const n of [1, 4, 7]) {
      await sk.click_atom(n, true, 'shift');
      await sk.click_bond(n, true, 'shift');
    }
    await sk.click_button('explicit_h');
    await checkpoint(page, 'add_explicit_hydrogens');
    await sk.click_button('explicit_h');
    await checkpoint(page, 'remove_explicit_hydrogens');
    await sk.click_tool('up');
    await checkpoint(page, 'change_bonds_to_up');
    await sk.click_tool('single');
    await checkpoint(page, 'change_bonds_to_single');
    await sk.click_button('plus_charge');
    await checkpoint(page, 'increase_charge');
    await sk.click_button('minus_charge');
    await checkpoint(page, 'decrease_charge');
    for (const tool of ['A', 'AH', 'Q', 'QH', 'M', 'MH', 'X', 'XH']) {
      await sk.click_tool(tool, true);
      await checkpoint(page, `change_atoms_to_${tool}`);
    }
    await sk.click_tool('down');
    await checkpoint(page, 'change_bonds_to_down');
    await sk.click_tool('single_either', true);
    await checkpoint(page, 'change_bonds_to_single_either');
    await sk.click_tool('double_either', true);
    await checkpoint(page, 'change_bonds_to_double_either');
    for (const tool of ['coordinate', 'double', 'triple', 'zero']) {
      await sk.click_tool(tool, true);
      await checkpoint(page, `change_bonds_to_${tool}`);
    }
    for (const tool of ['aromatic', 'any', 'single_double', 'single_aromatic', 'double_aromatic']) {
      await sk.click_tool(tool, true);
      await checkpoint(page, `change_bonds_to_${tool}`);
    }
    await sk.click_tool('Pd', true);
    await checkpoint(page, 'change_atoms_to_Pd_again');
    await sk.click_tool('C');
    await checkpoint(page, 'change_back_to_C');
    await sk.click_tool('last_picked_element');
    await checkpoint(page, 'last_picked_element');

    // Edit or Replace selected structure via the real right-click menu.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    for (const n of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    const copyText = await sk.selection_context_copy_text({ type: 'atom', index: 1 });
    expect(copyText).not.toBe('');
    for (const format of ['sdf', 'smi', 'cxsmi', 'inchi', 'inchikey', 'pdb']) {
      const copied = await sk.selection_context_copy_as_text({ type: 'atom', index: 1 }, format);
      expect(copied).not.toBe('');
    }
    await sk.selection_context_menu({ type: 'atom', index: 1 }, 'flip');
    await checkpoint(page, 'flip_context');
    await sk.click_button('undo');
    await sk.selection_context_menu({ type: 'atom', index: 1 }, 'cut');
    await checkpoint(page, 'cut_context');

    // Source-ordered Modify Atoms context-menu actions.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    for (const n of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    const atomTarget = { type: 'atom', index: 1 };
    for (const [action, checkpointName] of [
      ['add_explicit_hydrogens', 'add_explicit_hydrogens_context'],
      ['remove_explicit_hydrogens', 'remove_explicit_hydrogens_context'],
      ['+_charge', 'increase_charge_context'],
      ['–_charge', 'decrease_charge_context'],
      ['add_unpaired_electron', 'add_unpaired_electron_context'],
      ['remove_unpaired_electron', 'remove_unpaired_electron_context'],
    ]) {
      await sk.selection_context_menu(atomTarget, 'modify_atoms', action);
      await checkpoint(page, checkpointName);
    }
    for (const tool of ['A', 'AH', 'Q', 'QH', 'M', 'MH', 'X', 'XH']) {
      await sk.selection_context_menu(atomTarget, 'modify_atoms', 'replace_atoms_with', 'wildcard', tool);
      await checkpoint(page, `change_atoms_to_${tool}_context`);
    }
    // R-group placement has no deterministic Squish image reference
    // (SKETCH-1905), but retain the same real context-menu interactions.
    await sk.click_button('clear_selection');
    for (const n of [3, 4, 5, 6, 7, 8]) {
      await sk.click_atom(n, true, 'shift');
      await sk.click_bond(n, true, 'shift');
    }
    await sk.selection_context_menu({ type: 'atom', index: 3 }, 'modify_atoms', 'replace_atoms_with', 'new_r-group');
    await sk.click_button('clear_selection');
    await sk.click_tool('rect_btn');
    for (const n of [1, 2]) {
      await sk.click_atom(n, true, 'shift');
      await sk.click_bond(n, true, 'shift');
    }
    await sk.selection_context_menu(atomTarget, 'modify_atoms', 'replace_atoms_with', 'existing_r-group', 'r1');
    await sk.click_button('clear_selection');
    for (const n of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    await sk.open_selection_context_submenu(atomTarget, 'modify_atoms', 'set_element');
    await sendWidgetMousePress(page, 'periodic_table_btn');
    await page.waitForTimeout(25);
    await clickWidget(page, 'pd_btn');
    await closeActiveQtPopups(page);

    // Source-ordered Modify Bonds context-menu actions.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    for (const n of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    const bondTarget = { type: 'atom', index: 1 };
    for (const bond of ['single', 'double', 'triple', 'aromatic', 'up', 'down']) {
      await sk.selection_context_menu(bondTarget, 'modify_bonds', bond);
      await checkpoint(page, `${bond}_context`);
    }
    for (const bond of ['coordinate', 'zero_order', 'single_up/down', 'double_cis/trans']) {
      await sk.selection_context_menu(bondTarget, 'modify_bonds', 'other_type', bond);
      await checkpoint(page, `${bond.replace('/', '_')}_context`);
    }
    for (const bond of ['any', 'single/double', 'single/aromatic', 'double/aromatic']) {
      await sk.selection_context_menu(bondTarget, 'modify_bonds', 'query', bond);
      await checkpoint(page, `${bond.replace('/', '_')}_context`);
    }
  });
});
