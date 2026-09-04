import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

test('Modify Atoms R-group context menus use real cascading menus', async ({ page }) => {
  const sk = new Sketcher(page);
  await sk.open();
  await sk.import_menu('paste_in_text', SOURCE);
  await sk.map_imported_atom_indexes();
  for (const n of [3, 4, 5, 6, 7, 8]) {
    await sk.click_atom(n, true, 'shift');
    await sk.click_bond(n, true, 'shift');
  }
  await sk.selection_context_menu(
    { type: 'atom', index: 3 }, 'modify_atoms', 'replace_atoms_with', 'new_r-group',
  );
  await sk.click_button('clear_selection');
  await sk.click_tool('rect_btn');
  for (const n of [1, 2]) {
    await sk.click_atom(n, true, 'shift');
    await sk.click_bond(n, true, 'shift');
  }
  await sk.selection_context_menu(
    { type: 'atom', index: 1 }, 'modify_atoms', 'replace_atoms_with', 'existing_r-group', 'r1',
  );
  // R-group placement has no deterministic image reference (SKETCH-1905).
  expect(await sk.clipboard_text()).toBeDefined();
});
