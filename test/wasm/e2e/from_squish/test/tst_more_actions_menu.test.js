import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function requireRenderedGeometryBridge(page) {
  const available = await page.evaluate(
    () =>
      typeof Module._sketcher_get_atom_rect === 'function' &&
      typeof Module._sketcher_get_bond_rect === 'function',
  );
  test.skip(!available, 'requires the WASM artifact with rendered-geometry test bridge support');
}

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

async function selectSourceAtomsAndBonds(sk) {
  // Exact source sequence from suite_molviewer/tst_more_actions_menu/test.py.
  for (const n of [1, 2, 3, 4, 5, 6]) {
    await sk.click_bond(n, true, 'shift');
    await sk.click_atom(n, true, 'shift');
  }
}

async function setupSourceSelection(sk) {
  await sk.load_structure_for_test(SOURCE_STRUCTURE);
  await selectSourceAtomsAndBonds(sk);
}

test.describe('tst_more_actions_menu', () => {
  test.beforeEach(async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await requireRenderedGeometryBridge(page);

    // Source setup is Import -> Paste in Text, then Clear -> build_structure.
    // The standalone Paste in Text dialog is an explicit WASM limitation, so
    // this fixture loader is setup only and not Import-menu coverage.
    await sk.load_structure_for_test(SOURCE_STRUCTURE);
  });

  test('add_hydrogens_ignores_selection', async ({ page }) => {
    const sk = new Sketcher(page);
    await selectSourceAtomsAndBonds(sk);
    await sk.more_actions_menu('modify_all', 'add_explicit_hydrogens');
    await checkpoint(page, 'add_hydrogens_ignores_selection');
    await sk.click_button('undo');
  });

  test('Ctrl_X_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await selectSourceAtomsAndBonds(sk);
    await sk.type_text('sketcher_area', '<Ctrl+X>');
    await checkpoint(page, 'Ctrl_X_hotkey');
  });

  test('Ctrl_V_hotkey', async () => {
    test.skip(
      true,
      'Qt/WASM does not receive browser-delivered Ctrl+V; menu Paste remains covered separately.',
    );
  });

  test('Ctrl_A_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    // This is the source state immediately after Ctrl+V.  Browser Ctrl+V is
    // unavailable to Qt/WASM, so fixture setup replaces only that transition.
    await sk.load_structure_for_test(SOURCE_STRUCTURE);
    await sk.type_text('sketcher_area', '<Ctrl+A>');
    await checkpoint(page, 'Ctrl_A_hotkey');
  });

  test('Ctrl_Z_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);
    await sk.type_text('sketcher_area', '<Ctrl+X>');
    await sk.type_text('sketcher_area', '<Ctrl+Z>');
    await checkpoint(page, 'Ctrl_Z_hotkey');
  });

  test('Ctrl_Shift_Z_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);
    await sk.type_text('sketcher_area', '<Ctrl+X>');
    await sk.type_text('sketcher_area', '<Ctrl+Z>');
    await sk.type_text('sketcher_area', '<Ctrl+Shift+Z>');
    await checkpoint(page, 'Ctrl_Shift_Z_hotkey');
  });

  test('Ctrl_D_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);
    await sk.type_text('sketcher_area', '<Ctrl+X>');
    await sk.type_text('sketcher_area', '<Ctrl+Z>');
    await sk.type_text('sketcher_area', '<Ctrl+Shift+Z>');
    await sk.type_text('sketcher_area', '<Ctrl+D>');
    await checkpoint(page, 'Ctrl_D_hotkey');
  });

  test('Ctrl_I_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);
    await sk.type_text('sketcher_area', '<Ctrl+I>');
    await checkpoint(page, 'Ctrl_I_hotkey');
  });

  test('Ctrl_C_hotkey', async () => {
    test.skip(true, 'Qt/WASM does not expose application clipboard contents to the browser.');
  });

  test('Ctrl_F_hotkey', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);
    await sk.click_button('clear_selection');
    await sk.mouse_drag(0, 0, 100, 100, null, 'right');
    await sk.type_text('sketcher_area', '<Ctrl+F>');
    await checkpoint(page, 'Ctrl_F_hotkey');
  });

  test('more_actions_menu', async ({ page }) => {
    const sk = new Sketcher(page);
    await setupSourceSelection(sk);

    for (const action of [
      'flip_horizontal',
      'flip_vertical',
      'add_explicit_hydrogens',
      'remove_explicit_hydrogens',
    ]) {
      await sk.more_actions_menu('modify_all', action);
      await checkpoint(page, action);
    }

    // In the standalone build, Modify All retains the prior selection and
    // presents "Copy As".  The source then invokes "Copy All As", so clear
    // that retained selection with the real UI before its whole-structure
    // export block.
    await sk.click_button('clear_selection');
    for (const format of ['sdf', 'smi', 'cxsmi', 'inchi', 'inchikey', 'pdb']) {
      await sk.more_actions_menu('copy_all_as', format);
      // Squish records getClipboardText() here.  Qt/WASM performs the real
      // menu click, but its QApplication clipboard is not browser-readable;
      // see the explicit skipped clipboard-reference case below.
    }

    // Keep this stateful order identical to the source buttons_list.
    for (const [action, reference] of [
      ['fit_to_screen', 'fit_to_screen'],
      ['copy_all', null],
      ['select_all', 'select_all_1'],
      ['clear_selection', 'clear_selection'],
      ['select_all', 'select_all_2'],
      ['invert_selection', 'invert_selection'],
      ['select_all', 'select_all_3'],
      ['cut', null],
      ['paste', 'paste'],
      ['undo', null],
      ['redo', 'redo'],
    ]) {
      await sk.more_actions_menu(action);
      if (reference === null) {
        // These are clipboard-only Squish references.  The menu actions are
        // still exercised with real mouse input, but text verification is
        // unavailable in Qt/WASM.
      } else {
        await checkpoint(page, reference);
      }
    }
  });
});
