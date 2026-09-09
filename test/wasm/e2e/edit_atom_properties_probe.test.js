import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
test('edit atom properties dialog inventory', async ({ page }) => {
  const sk = new Sketcher(page); await sk.open();
  await sk.import_menu('paste_in_text', 'CCC[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)CC');
  await sk.map_imported_atom_indexes(); await sk.click_atom(4, true, 'shift');
  await sk.selection_context_menu({ type: 'atom', index: 4 }, 'edit_atom_properties');
  // Replace the temporary widget inventory throw with stable assertions on
  // the visible Atom dialog defaults, then dismiss it as a user would.
  expect((await sk.widget_state('atom_element_le')).text).toBe('C');
  expect((await sk.widget_state('atom_isotope_sb')).text).toBe('');
  expect((await sk.widget_state('atom_charge_sb')).text).toBe('0');
  await sk.edit_atom_properties({ click_ok: false, click_cancel: true });
});
