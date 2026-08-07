import { expect, test } from '@playwright/test';
import {
  clickMenuAction,
  clickWidget,
  clipboardText,
  exportText,
  focusCanvas,
  importText,
  openSketcher,
  setClipboardText,
} from '../wrappers/sketcher_wasm.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function requireTestBridge(page) {
  const available = await page.evaluate(
    () => typeof Module._sketcher_get_menu_action_rect === 'function',
  );
  test.skip(!available, 'requires a WASM artifact rebuilt with playwright_test_bridge.cpp');
}

async function moreAction(page, action) {
  // Opening More Actions first preserves production menu-state behavior.
  await clickWidget(page, 'more_actions_btn');
  if (
    [
      'Flip Horizontal',
      'Flip Vertical',
      'Add Explicit Hydrogens',
      'Remove Explicit Hydrogens',
    ].includes(action)
  ) {
    await clickMenuAction(page, 'Modify All');
  }
  await clickMenuAction(page, action);
}

async function expectCanvasCheckpoint(page, name) {
  await page.mouse.move(0, 0);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.beforeEach(async ({ page }) => {
  await openSketcher(page);
  await requireTestBridge(page);
  await importText(page, SOURCE_STRUCTURE);
});

test.describe('ported tst_more_actions_menu', () => {
  test('modify-all actions preserve valid editable structures', async ({ page }) => {
    // Squish: modify_all -> flip_horizontal, flip_vertical, add/remove H.
    for (const [action, checkpoint] of [
      ['Flip Horizontal', 'more-actions-flip-horizontal'],
      ['Flip Vertical', 'more-actions-flip-vertical'],
      ['Add Explicit Hydrogens', 'more-actions-add-explicit-hydrogens'],
      ['Remove Explicit Hydrogens', 'more-actions-remove-explicit-hydrogens'],
    ]) {
      await moreAction(page, action);
      await expect.poll(() => exportText(page)).not.toBe('');
      await expectCanvasCheckpoint(page, checkpoint);
    }
  });

  test('keyboard cut, paste, selection, undo, redo, and fit match menu actions', async ({
    page,
  }) => {
    // Squish checkpoint sequence: Ctrl+X, Ctrl+V, Ctrl+A, Ctrl+Z,
    // Ctrl+Shift+Z, Ctrl+D, Ctrl+I, Ctrl+C, Ctrl+F.
    await focusCanvas(page);
    // The Squish test first selects atoms and bonds before cutting. Selecting
    // the whole structure gives the standalone web build equivalent state.
    await page.keyboard.press('ControlOrMeta+a');
    await page.keyboard.press('ControlOrMeta+x');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-x');

    await page.keyboard.press('ControlOrMeta+v');
    await expect.poll(() => exportText(page)).not.toBe('');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-v');

    await page.keyboard.press('ControlOrMeta+a');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-a');

    await page.keyboard.press('ControlOrMeta+z');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-z');

    await page.keyboard.press('ControlOrMeta+Shift+z');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-shift-z');

    await page.keyboard.press('ControlOrMeta+d');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-d');

    await page.keyboard.press('ControlOrMeta+i');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-i');

    await page.keyboard.press('ControlOrMeta+a');
    await page.keyboard.press('ControlOrMeta+c');
    expect(await clipboardText(page)).not.toBe('');

    // Squish clears selection and moves the structure before fitting it.
    await page.keyboard.press('Escape');
    const canvas = page.locator('#screen canvas');
    const box = await canvas.boundingBox();
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    await page.mouse.down({ button: 'right' });
    await page.mouse.move(box.x + box.width / 2 + 100, box.y + box.height / 2 + 100);
    await page.mouse.up({ button: 'right' });
    await page.keyboard.press('ControlOrMeta+f');
    await expectCanvasCheckpoint(page, 'more-actions-ctrl-f');
  });

  test('More Actions menu commands mirror the original command sequence', async ({ page }) => {
    // The source test invokes these entries in this order.  Keep the same
    // order because cut/copy/selection state controls action availability.
    for (const [action, checkpoint] of [
      ['Fit to Screen', 'more-actions-fit-to-screen'],
      ['Select All', 'more-actions-select-all-1'],
      ['Clear Selection', 'more-actions-clear-selection'],
      ['Select All', 'more-actions-select-all-2'],
      ['Invert Selection', 'more-actions-invert-selection'],
      ['Select All', 'more-actions-select-all-3'],
      ['Cut', 'more-actions-cut'],
      ['Paste', 'more-actions-paste'],
      ['Undo', 'more-actions-undo'],
      ['Redo', 'more-actions-redo'],
    ]) {
      if (action === 'Paste') {
        await setClipboardText(page, SOURCE_STRUCTURE);
      }
      await moreAction(page, action);
      await expectCanvasCheckpoint(page, checkpoint);
    }
  });
});
