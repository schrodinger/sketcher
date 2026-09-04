import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

test('Modify Atoms context actions use real cascading menus', async ({ page }) => {
  const sk = new Sketcher(page);
  await sk.open();
  await sk.import_menu('paste_in_text', SOURCE);
  await sk.map_imported_atom_indexes();
  for (const n of [1, 2, 3, 4, 5, 6]) {
    await sk.click_bond(n, true, 'shift');
    await sk.click_atom(n, true, 'shift');
  }
  const target = { type: 'atom', index: 1 };
  await sk.selection_context_menu(target, 'modify_atoms', 'add_explicit_hydrogens');
  await sk.selection_context_menu(target, 'modify_atoms', 'remove_explicit_hydrogens');
  await sk.selection_context_menu(target, 'modify_atoms', 'replace_atoms_with', 'wildcard', 'XH');
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot('modify-atoms-actions.png');
});
