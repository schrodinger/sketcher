import { expect, test } from '@playwright/test';
import {
  clickWidget,
  drawingAreaCenter,
  drawBond,
  exportText,
  focusCanvas,
  isEmpty,
  openSketcher,
} from '../wrappers/sketcher_wasm.js';

test.beforeEach(async ({ page }) => {
  await openSketcher(page);
});

test.describe('ported Squish canvas editing', () => {
  // Source coverage: tst_erase_mode, tst_move_mode, and tst_select_mode_*.
  test('erase tool removes a placed atom', async ({ page }) => {
    const center = await drawingAreaCenter(page);
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => isEmpty(page)).toBe(false);

    await clickWidget(page, 'erase_btn');
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => isEmpty(page)).toBe(true);
  });

  test('undo and redo restore a drawn bond', async ({ page }) => {
    await drawBond(page);
    await expect.poll(() => exportText(page)).toBe('CC');

    await focusCanvas(page);
    await page.keyboard.press('ControlOrMeta+z');
    await expect.poll(() => isEmpty(page)).toBe(true);

    await page.keyboard.press('ControlOrMeta+Shift+z');
    await expect.poll(() => exportText(page)).toBe('CC');
  });

  test('increase-charge tool changes the selected atom', async ({ page }) => {
    const center = await drawingAreaCenter(page);
    await page.mouse.click(center.x, center.y);
    await clickWidget(page, 'increase_charge_btn');
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => exportText(page)).toContain('+');
  });
});
