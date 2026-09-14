import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getExportedHelm,
  selectAll,
  clickWidget,
  requireRect,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Monomer Mutation', () => {
  test('mutate selected peptide monomers to cysteine', async ({ page }) => {
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

  test('incompatible monomer button is disabled for a peptide selection', async ({ page }) => {
    await clickWidget(page, 'monomeric_btn');
    await clickWidget(page, 'amino_monomer_btn');

    await page.evaluate(() => {
      Module.sketcher_import_text('PEPTIDE1{A.G.L}$$$$V2.0');
    });
    await clickWidget(page, 'fit_btn');

    await selectAll(page);

    // Switch to nucleic acid page — peptide selection should persist
    await clickWidget(page, 'nucleic_monomer_btn');

    // A nucleic acid base can't be applied to a peptide selection, so the
    // button must be disabled. Assert that directly rather than clicking it and
    // checking the molecule is unchanged, which would also pass if the button
    // were enabled but the mutation silently did nothing.
    await expect
      .poll(async () => (await requireRect(page, 'widget:na_u_btn')).enabled, { timeout: 5000 })
      .toBe(false);

    const helm = await getExportedHelm(page);
    expect(helm).toContain('PEPTIDE1{A.G.L}');
  });

  test('undo restores original monomers after mutation', async ({ page }) => {
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
