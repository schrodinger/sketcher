import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
const SOURCE = 'CCC[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)CC';
test('tst_edit_atom_properties atom tab', async ({ page }) => {
  const sk = new Sketcher(page); await sk.open();
  await sk.import_menu('paste_in_text', SOURCE); await sk.map_imported_atom_indexes();
  await sk.click_atom(4, true, 'shift');
  await sk.selection_context_menu({ type: 'atom', index: 4 }, 'edit_atom_properties');
  expect((await sk.widget_state('element_le')).text).toBe('C');
  await sk.edit_atom_properties({ element: 'O' });
  await expect(page.locator('#screen canvas')).toHaveScreenshot('Change_element_to_O.png');
});
