import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { hideMouseMarker } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
const SGROUP_ATOMS = [4, 5, 6, 7];
const SGROUP_XBONDS = [4, 5, 6];

test.setTimeout(180_000);

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_bracket_subgroup', () => {
  // Temporary focused regression check for the 2026-4 COP bracket report.
  // It deliberately uses the same real pointer/keyboard wrapper sequence as
  // the Squish port and the established pre-regression COP/HT visual baseline.
  test('temporary COP brackets regression', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    await sk.add_sgroup('COP', SGROUP_ATOMS, SGROUP_XBONDS, 'HT', '');
    await page.mouse.move(0, 0);
    await hideMouseMarker(page);
    await expect(page.locator('#screen canvas')).toHaveScreenshot('bracket-subgroup-COP-HT.png');
  });

  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source imports, captures its structure coordinates, clears, and rebuilds.
    // Retain the imported structure and map source indexes to its live browser
    // geometry so the same selection/menu workflow is not desktop-pixel bound.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();

    for (const [type, connect, label, name] of [
      ['SRU', 'HT', '', 'bracket_subgroup_SRU_HT_nolabel'],
      ['SRU', 'HH', '', 'bracket_subgroup_SRU_HH_nolabel'],
      ['SRU', 'EU', '', 'bracket_subgroup_SRU_EU_nolabel'],
      ['SRU', 'HT', '4', 'bracket_subgroup_SRU_HT_constantlabel'],
      ['SRU', 'HT', '3-7', 'bracket_subgroup_SRU_HT_rangelabel'],
      ['COP', 'HT', '', 'bracket_subgroup_COP_HT'],
      ['COP', 'HH', '', 'bracket_subgroup_COP_HH'],
      ['COP', 'EU', '', 'bracket_subgroup_COP_EU'],
    ]) {
      await sk.add_sgroup(type, SGROUP_ATOMS, SGROUP_XBONDS, connect, label);
      await checkpoint(page, name);
      // The source undoes every case except its final COP/EU reference.
      if (name !== 'bracket_subgroup_COP_EU') await sk.click_button('undo');
    }
  });
});
