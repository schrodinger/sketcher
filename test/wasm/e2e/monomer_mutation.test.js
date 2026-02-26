import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getExportedHelm,
  enableMonomericMode,
  selectAll,
  clickWidget,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Monomer Mutation', () => {
  test('mutate selected peptide monomers to cysteine', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{A.G.L}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    await selectAll(page);
    await clickWidget(page, 'cys_btn');

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('PEPTIDE1{C.C.C}');
  });

  test('mutate selected nucleic acid bases only', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'nucleic_monomer_btn');

    await page.evaluate(() => {
      Module.sketcher_import_text('RNA1{R(A)P.R(G)P.R(C)P}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    await selectAll(page);
    await clickWidget(page, 'na_u_btn');

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('R(U)P');

    const helm = await getExportedHelm(page);
    expect(helm).not.toContain('R(A)');
    expect(helm).not.toContain('R(G)');
    expect(helm).not.toContain('R(C)');
  });

  test('clicking disabled incompatible monomer button does not mutate', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{A.G.L}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    await selectAll(page);

    // Switch to nucleic acid page — peptide selection should persist
    await clickWidget(page, 'nucleic_monomer_btn');

    // Click a nucleic acid base button (should be disabled for peptide selection)
    await clickWidget(page, 'na_u_btn');

    // Verify no mutation occurred
    const helm = await getExportedHelm(page);
    expect(helm).toContain('PEPTIDE1{A.G.L}');
  });

  test('undo restores original monomers after mutation', async ({ page }) => {
    await enableMonomericMode(page);
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{A.G.L}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    await selectAll(page);
    await clickWidget(page, 'cys_btn');

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('PEPTIDE1{C.C.C}');

    // Undo the mutation
    await focusCanvas(page);
    await page.keyboard.press('ControlOrMeta+z');

    await expect.poll(() => getExportedHelm(page), { timeout: 5000 }).toContain('PEPTIDE1{A.G.L}');
  });
});
