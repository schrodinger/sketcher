import { expect, test } from '@playwright/test';
import {
  clickMenuAction,
  clickWidget,
  drawingAreaCenter,
  exportText,
  loadStructureForTest,
  mouseDrag,
  openSketcher,
} from '../wrappers/sketcher_wasm.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function requireTestBridge(page) {
  const available = await page.evaluate(
    () => typeof Module._sketcher_get_menu_action_rect === 'function',
  );
  test.skip(!available, 'requires a WASM artifact rebuilt with playwright_test_bridge.cpp');
}

async function moreAction(page, action) {
  await clickWidget(page, 'more_actions_btn');
  await clickMenuAction(page, action);
}

test.describe('ported tst_move_mode and tst_select_mode_active_selection', () => {
  test.beforeEach(async ({ page }) => {
    await openSketcher(page);
    await requireTestBridge(page);
    await loadStructureForTest(page, SOURCE_STRUCTURE);
  });

  test('Move/Rotate preserves the selected structure during a pointer drag', async ({ page }) => {
    // Squish selects atoms/bonds, changes to Move/Rotate, then drags a selected
    // atom. Selecting all is the standalone equivalent of that active state.
    await moreAction(page, 'Select All');
    await clickWidget(page, 'move_rotate_btn');

    const center = await drawingAreaCenter(page);
    await mouseDrag(page, center, { x: center.x + 70, y: center.y + 55 });

    await expect.poll(() => exportText(page)).not.toBe('');
    await page.mouse.move(0, 0);
    await expect(page.locator('#screen canvas')).toHaveScreenshot('move-selected-structure.png');
  });
});
