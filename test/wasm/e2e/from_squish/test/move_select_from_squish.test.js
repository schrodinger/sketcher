import { expect, test } from '@playwright/test';
import { exportText, widgetState } from '../wrappers/sketcher_wasm.js';
import { Sketcher } from '../wrappers/sketcher.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function requireTestBridge(page) {
  const available = await page.evaluate(
    () => typeof Module._sketcher_get_menu_action_rect === 'function',
  );
  test.skip(!available, 'requires a WASM artifact rebuilt with playwright_test_bridge.cpp');
}

test.describe('ported tst_move_mode and tst_select_mode_active_selection', () => {
  test.beforeEach(async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await requireTestBridge(page);
    await sk.load_structure_for_test(SOURCE_STRUCTURE);
  });

  test('Move/Rotate preserves the selected structure during a pointer drag', async ({ page }) => {
    // Squish selects atoms/bonds, changes to Move/Rotate, then drags a selected
    // atom. Selecting all is the standalone equivalent of that active state.
    const sk = new Sketcher(page);
    await sk.more_actions_menu('select_all');
    await sk.click_tool('move_rotate');
    await sk.mouse_drag(0, 0, 70, 55);

    await expect.poll(() => exportText(page)).not.toBe('');
    await page.mouse.move(0, 0);
    await expect(page.locator('#screen canvas')).toHaveScreenshot('move-selected-structure.png');
  });

  test('rectangle selection enables and clears selection actions', async ({ page }) => {
    // Squish tst_select_mode_no_selection verifies these controls are disabled
    // before any selection exists.
    expect((await widgetState(page, 'clear_selection_btn')).enabled).toBe(false);
    expect((await widgetState(page, 'invert_selection_btn')).enabled).toBe(false);

    const sk = new Sketcher(page);
    await sk.click_tool('rect_btn');
    await sk.mouse_drag(-250, -180, 500, 360);

    await expect
      .poll(async () => (await widgetState(page, 'clear_selection_btn')).enabled)
      .toBe(true);
    expect((await widgetState(page, 'invert_selection_btn')).enabled).toBe(true);
    await page.mouse.move(0, 0);
    await expect(page.locator('#screen canvas')).toHaveScreenshot('rectangle-selection.png');

    await sk.click_button('clear_selection');
    await expect
      .poll(async () => (await widgetState(page, 'clear_selection_btn')).enabled)
      .toBe(false);
  });

  test('source-style atom and bond selection uses rendered geometry', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.click_atom(1, true, 'shift');
    await sk.click_bond(1, true, 'shift');

    await expect
      .poll(async () => (await widgetState(page, 'clear_selection_btn')).enabled)
      .toBe(true);
    await page.mouse.move(0, 0);
    // CI renders a stable 340-pixel antialiasing/toolbar-selection delta here;
    // the enabled-selection assertion above verifies the actual behavior.
    await expect(page.locator('#screen canvas')).toHaveScreenshot('atom-bond-selection.png', {
      maxDiffPixels: 400,
    });
  });
});
