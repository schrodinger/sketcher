import { expect, test } from '@playwright/test';
import {
  clickPopupTool,
  contextMenuAction,
  getExportedSmiles,
  loadStructure,
  waitForSketcherReady,
  widgetState,
} from './e2e_helpers.js';

// SKETCH-1652: on WASM, closing a Qt::Popup window (a context menu or a side
// bar tool popup) used to leave keyboard focus somewhere other than the
// sketcher, so shortcuts were ignored until the user clicked the scene. Each
// test presses its shortcut straight after the popup closes, as a user would.
test.describe('Keyboard shortcuts after a popup closes', () => {
  test.beforeEach(async ({ page }) => {
    await waitForSketcherReady(page);
    await loadStructure(page, 'CC');
  });

  test('undo works right after changing a bond type from the context menu', async ({ page }) => {
    await contextMenuAction(page, 'bond:0', 'Double');
    await expect.poll(() => getExportedSmiles(page)).toBe('C=C');

    await page.keyboard.press('ControlOrMeta+z');
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('CC');
  });

  test('select all works right after picking a side bar tool from its popup', async ({ page }) => {
    await clickPopupTool(page, 'triple_btn');

    // The fit button is restyled whenever there is a selection
    await page.keyboard.press('ControlOrMeta+a');
    await expect
      .poll(async () => (await widgetState(page, 'fit_btn')).styleSheet, { timeout: 5000 })
      .not.toBe('');
  });
});
