import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getCanvasCenter,
  getExportedSmiles,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Drawing', () => {
  test('click on empty canvas places a carbon atom', async ({ page }) => {
    const center = await getCanvasCenter(page);

    await page.mouse.click(center.x, center.y);

    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('C');
  });

  // Each element shortcut gets its own test with a fresh page to avoid
  // interference from programmatic clears, which can leave the sketcher
  // in a transient state that ignores subsequent mouse clicks.
  for (const { key, smiles: expected } of [
    { key: 'N', smiles: 'N' },
    { key: 'O', smiles: 'O' },
    { key: 'S', smiles: 'S' },
  ]) {
    test(`keyboard shortcut "${key}" places ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await focusCanvas(page);
      await page.keyboard.press(key);
      await page.mouse.click(center.x, center.y);

      await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe(expected);
    });
  }

  test('click-drag draws a bond', async ({ page }) => {
    const center = await getCanvasCenter(page);
    const dragDistance = 100;

    await page.mouse.move(center.x, center.y);
    await page.mouse.down();
    await page.mouse.move(center.x + dragDistance, center.y, { steps: 10 });
    await page.mouse.up();

    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('CC');
  });

  test('draw bond screenshot', async ({ page }) => {
    const center = await getCanvasCenter(page);

    await page.mouse.move(center.x, center.y);
    await page.mouse.down();
    await page.mouse.move(center.x + 100, center.y, { steps: 10 });
    await page.mouse.up();

    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('CC');

    // Move the mouse away so the cursor doesn't affect the screenshot
    await page.mouse.move(0, 0);

    const canvas = page.locator('#screen canvas');
    await expect(canvas).toHaveScreenshot('e2e-draw-bond.png');
  });
});
