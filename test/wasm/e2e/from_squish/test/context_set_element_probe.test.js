import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { clickWidget, closeActiveQtPopups, sendWidgetMousePress } from '../wrappers/sketcher_wasm.js';

test('Modify Atoms Set Element opens the periodic-table popup by mouse', async ({ page }) => {
  const sk = new Sketcher(page);
  await sk.open();
  await sk.import_menu('paste_in_text', 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl');
  await sk.map_imported_atom_indexes();
  await sk.click_button('select_all');
  await sk.open_selection_context_submenu({ type: 'atom', index: 1 }, 'modify_atoms', 'set_element');
  await sendWidgetMousePress(page, 'periodic_table_btn');
  await page.waitForTimeout(25);
  await clickWidget(page, 'pd_btn');
  await closeActiveQtPopups(page);
  await page.waitForTimeout(100);
  expect(await sk.copy_all_as_text('sdf')).toContain(' Pd ');
});
