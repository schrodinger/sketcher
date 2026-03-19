import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  getDrawingAreaCenter,
  getExportedSmiles,
  clickWidget,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Import', () => {
  test('import molecule via API then extend by drawing', async ({ page }) => {
    // Import a molecule programmatically (the LiveDesign workflow)
    await page.evaluate(() => {
      Module.sketcher_import_text('C');
    });
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('C');

    // Fit the view so the imported atom is centered in the drawing area
    await clickWidget(page, 'fit_btn');

    // Use the drawing area center (not canvas center) to hit the atom
    const center = await getDrawingAreaCenter(page);

    // Draw a bond from the atom outward to extend the molecule
    await clickWidget(page, 'single_bond_btn');
    await page.mouse.move(center.x, center.y);
    await page.mouse.down();
    await page.mouse.move(center.x + 100, center.y, { steps: 10 });
    await page.mouse.up();

    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('CC');
  });
});
