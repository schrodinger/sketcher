import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getCanvasCenter,
  getExportedSmiles,
  isSketcherEmpty,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Undo / Redo', () => {
  test('undo and redo', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place an atom
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('C');

    // Undo — sketcher should be empty again
    await focusCanvas(page);
    await page.keyboard.press('ControlOrMeta+z');
    await expect.poll(() => isSketcherEmpty(page), { timeout: 5000 }).toBe(true);

    // Redo — atom should reappear
    await focusCanvas(page);
    await page.keyboard.press('ControlOrMeta+Shift+z');
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('C');
  });
});
