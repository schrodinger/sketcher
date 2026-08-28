import { expect, test } from '@playwright/test';
import {
  clearSketcher,
  exportText,
  importText,
  isEmpty,
  openSketcher,
} from '../wrappers/sketcher_wasm.js';

test.beforeEach(async ({ page }) => {
  await openSketcher(page);
});

test.describe.skip('ported Squish import and export', () => {
  // Re-enable after the standalone WASM import/export dialogs have stable
  // browser-side geometry for real click-and-type interaction. These tests
  // must not use the fixture loader because they cover the UI flow itself.
  // Source coverage: tst_import_menu, tst_export_menu, and tst_miscellaneous.
  test('imports SMILES and exports the current structure', async ({ page }) => {
    await importText(page, 'C=O');
    await expect.poll(() => exportText(page)).toBe('C=O');
  });

  test('imports a molfile and exports SMILES', async ({ page }) => {
    const molfile = `\n  Sketcher          2D\n\n  2  1  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\nM  END`;
    await importText(page, molfile);
    await expect.poll(() => exportText(page)).toBe('C=O');
  });

  test('clear removes imported content', async ({ page }) => {
    await importText(page, 'c1ccccc1');
    await expect.poll(() => isEmpty(page)).toBe(false);
    await clearSketcher(page);
    await expect.poll(() => isEmpty(page)).toBe(true);
  });

  test('exports SVG image data', async ({ page }) => {
    await importText(page, 'C=O');
    const svg = await page.evaluate(() =>
      atob(Module.sketcher_export_image(Module.ImageFormat.SVG)),
    );
    expect(svg).toContain('<svg');
  });
});
