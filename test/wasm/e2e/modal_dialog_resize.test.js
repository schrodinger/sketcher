import { expect, test } from '@playwright/test';
import {
  activateAction,
  closeDialog,
  expectDialog,
  getRect,
  waitForSketcherReady,
} from './e2e_helpers.js';

const TITLE = 'Paste in Structure Text';

/**
 * Return the {x, y, width, height} of the foreground modal dialog. The bridge
 * reports it relative to the sketcher, whose top-left is the page's top-left.
 */
async function dialogGeometry(page) {
  const rect = await getRect(page, 'dialog:title');
  return rect && { x: rect.x, y: rect.y, width: rect.width, height: rect.height };
}

/**
 * Resize the page, and with it the canvas, then wait for the dialog to reach
 * the expected position. Qt only sees the new canvas size once the browser's
 * ResizeObserver fires, so the move happens asynchronously.
 */
async function resizeAndExpectDialogAt(page, width, height, size, x, y) {
  await page.setViewportSize({ width, height });
  await expect.poll(() => dialogGeometry(page)).toEqual({ x, y, ...size });
}

test.describe('modal dialogs on canvas resize', () => {
  test('keeps an open modal dialog on the canvas as it shrinks', async ({ page }) => {
    await waitForSketcherReady(page);
    await activateAction(page, 'Paste in Text...');
    await expectDialog(page, TITLE);

    const initial = await dialogGeometry(page);
    const size = { width: initial.width, height: initial.height };
    const viewport = page.viewportSize();
    // The dialog starts centered on the larger canvas
    expect(initial.x + initial.width).toBeLessThan(viewport.width);
    expect(initial.y + initial.height).toBeLessThan(viewport.height);

    // Shrink the canvas so the centered dialog sticks out past its bottom-right
    // corner. The dialog moves just far enough to fit.
    const margin = 40;
    const fitWidth = initial.width + margin;
    const fitHeight = initial.height + margin;
    expect(initial.x + initial.width).toBeGreaterThan(fitWidth);
    expect(initial.y + initial.height).toBeGreaterThan(fitHeight);
    await resizeAndExpectDialogAt(page, fitWidth, fitHeight, size, margin, margin);

    // Growing the canvas leaves a dialog that already fits where it is
    await page.setViewportSize(viewport);
    await expect
      .poll(() => getRect(page, 'widget:view').then((r) => r?.width))
      .toBeGreaterThan(fitWidth);
    expect(await dialogGeometry(page)).toEqual({ x: margin, y: margin, ...size });

    // A canvas smaller than the dialog keeps the dialog's title bar on screen
    await resizeAndExpectDialogAt(
      page,
      Math.round(initial.width / 2),
      Math.round(initial.height / 2),
      size,
      0,
      0,
    );

    // Once there's room again, the dialog is still modal and still usable
    await page.setViewportSize(viewport);
    await expect
      .poll(() => getRect(page, 'widget:view').then((r) => r?.width))
      .toBeGreaterThan(fitWidth);
    await closeDialog(page, 'Cancel');
    await expect.poll(() => getRect(page, 'dialog:title')).toBeNull();
  });
});
