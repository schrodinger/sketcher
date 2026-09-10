import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';

test.setTimeout(60_000);

test.describe('tst_query_bond_crash', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    let browserCrashed = false;
    page.on('crash', () => { browserCrashed = true; });
    await sk.open();

    // SKETCH-2399: repeatedly apply Single/Double to an imported bond.
    await sk.import_menu('paste_in_text', 'OC1CCCC[C@@H]1Cl');
    await sk.map_imported_atom_indexes();
    await sk.click_tool('single_double', true);
    for (let count = 1; count < 10; count += 1) {
      await sk.click_bond(1);
      await page.waitForTimeout(1_000);
      expect(browserCrashed).toBe(false);
    }

    // SKETCH-2409: applying Up twice through the selection context menu.
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl');
    await sk.map_imported_atom_indexes();
    await sk.click_bond(4, true);
    await sk.selection_context_menu({ type: 'bond', index: 4 }, 'up');
    await sk.selection_context_menu({ type: 'bond', index: 4 }, 'up');
    await page.waitForTimeout(5_000);
    expect(browserCrashed).toBe(false);
  });
});
