import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
test('wildcard context path', async ({ page }) => {
  const browserErrors = [];
  page.on('pageerror', (error) => browserErrors.push(`pageerror: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') browserErrors.push(`console: ${message.text()}`);
  });
  const sk = new Sketcher(page);
  try {
    await sk.open();
  } catch (error) {
    throw new Error(`Sketcher did not initialize: ${error.message}\n${browserErrors.join('\n')}`);
  }
  await sk.import_menu('paste_in_text', 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl');
  await sk.map_imported_atom_indexes();
  for (const n of [1,2,3,4,5,6]) await sk.click_atom(n, true, 'shift');
  // Use the same real right-click/cascading-menu path as the Squish wrapper.
  // The bridge inventory was useful while implementing it, but a deliberate
  // throw made this diagnostic probe a permanent suite failure.
  await sk.selection_context_menu(
    { type: 'atom', index: 1 }, 'modify_atoms', 'replace_atoms_with', 'wildcard', 'XH',
  );
  expect((await sk.widget_state('clear_selection_btn')).enabled).toBe(true);
  expect(browserErrors).toEqual([]);
});
