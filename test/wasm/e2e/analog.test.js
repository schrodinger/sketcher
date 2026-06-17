import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  getDrawingAreaCenter,
  getExportedHelm,
  enableMonomericMode,
  selectAll,
  clickWidget,
  clickPopupButton,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Amino Acid Analog Tests', () => {
  test('HELM round-trip preserves non-natural analog symbols', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');

    // Import HELM with non-natural analog symbols (brackets required for
    // multi-character HELM symbols)
    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{[dA].[meA].A}$$$$V2.0');
    });

    // Verify the analog symbols survive export
    const helm = await getExportedHelm(page);
    expect(helm).toContain('[dA]');
    expect(helm).toContain('[meA]');
  });

  test('select analog via popup and draw', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    // Click alanine to activate it (checks the button)
    await clickWidget(page, 'ala_btn');

    // Select the dA analog via the C++ API (popup buttons can't be
    // targeted by Playwright in WASM — see clickPopupButton helper).
    await clickPopupButton(page, 'analog_dA_btn');

    // Place the analog on canvas
    const center = await getDrawingAreaCenter(page);
    await page.mouse.click(center.x, center.y);

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('{[dA]}');
  });

  test('mutate non-natural analogs to standard amino acid', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    // Import a peptide chain with non-natural analogs
    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{[dA].[meA].[dC]}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    // Select all monomers
    await selectAll(page);

    // Click Glycine button to mutate all selected to G
    await clickWidget(page, 'gly_btn');

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('PEPTIDE1{G.G.G}');
  });
});
