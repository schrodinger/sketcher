import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
const SKIP_NESTED_MENUS = process.env.PLAYWRIGHT_SKIP_NESTED_MENUS === '1';

// This source-style workflow rebuilds the molecule twice with visible input.
// Preserve the fast default timeout while allowing paced replay inspection.
test.setTimeout(180_000);

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
  await hideMouseMarker(page);
  if (process.env.PLAYWRIGHT_SKIP_SCREENSHOTS === '1') {
    return;
  }
  const replaySuffix = process.env.PLAYWRIGHT_REBUILD_STRUCTURES === '1' ? '-replay' : '';
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}${replaySuffix}.png`);
}

async function selectSourceAtomsAndBonds(sk) {
  // Exact source sequence from suite_molviewer/tst_more_actions_menu/test.py.
  for (const n of [1, 2, 3, 4, 5, 6]) {
    await sk.click_bond(n, true, 'shift');
    await sk.click_atom(n, true, 'shift');
  }
}

test.describe('tst_more_actions_menu', () => {
  test('main', async ({ page }, testInfo) => {
    const sk = new Sketcher(page);
    await sk.open();
    await requireRenderedGeometryBridge(page);

    // Source setup is Import -> Paste in Text, then Clear -> build_structure.
    // The standalone Paste in Text dialog is covered by its own suite; fixture
    // setup here only establishes the same deterministic starting structure.
    await sk.load_structure_for_test(SOURCE_STRUCTURE);
    await selectSourceAtomsAndBonds(sk);

    if (!SKIP_NESTED_MENUS) {
      await test.step('add_hydrogens_ignores_selection', async () => {
        await sk.more_actions_menu('modify_all', 'add_explicit_hydrogens');
        await checkpoint(page, 'add_hydrogens_ignores_selection');
        await sk.click_button('undo');
      });
    } else {
      testInfo.annotations.push({
        type: 'temporary probe exclusion',
        description: 'Nested Modify All menu actions are disabled for diagnosis.',
      });
    }
    await test.step('Ctrl_X_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+X>');
      await checkpoint(page, 'Ctrl_X_hotkey');
    });

    // Browser-delivered Ctrl+V never reaches the Qt/WASM application. Clear
    // normally, then restore the exact post-paste source state as a fixture.
    await sk.click_button('clear');
    testInfo.annotations.push({
      type: 'WASM limitation',
      description: 'Ctrl_V_hotkey: browser Ctrl+V does not reach Qt/WASM.',
    });
    await sk.load_structure_for_test(SOURCE_STRUCTURE);

    await test.step('Ctrl_A_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+A>');
      await checkpoint(page, 'Ctrl_A_hotkey');
    });
    await test.step('Ctrl_Z_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+Z>');
      await checkpoint(page, 'Ctrl_Z_hotkey');
    });
    await test.step('Ctrl_Shift_Z_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+Shift+Z>');
      await checkpoint(page, 'Ctrl_Shift_Z_hotkey');
    });
    await test.step('Ctrl_D_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+D>');
      await checkpoint(page, 'Ctrl_D_hotkey');
    });

    // This is the source's second import -> rebuild -> selection setup.
    // `load_structure_for_test` imports into the current model, while the
    // source Import/Paste workflow replaces it before rebuilding.
    await sk.click_button('clear');
    await sk.load_structure_for_test(SOURCE_STRUCTURE);
    await selectSourceAtomsAndBonds(sk);

    await test.step('Ctrl_I_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+I>');
      await checkpoint(page, 'Ctrl_I_hotkey');
    });
    await test.step('Ctrl_C_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+C>');
      await expect(await sk.clipboard_text()).toContain('V3000');
    });
    await test.step('Ctrl_F_hotkey', async () => {
      await sk.click_button('clear_selection');
      await sk.mouse_drag(0, 0, 100, 100, null, 'right');
      await sk.type_text('sketcher_area', '<Ctrl+F>');
      await checkpoint(page, 'Ctrl_F_hotkey');
    });

    if (!SKIP_NESTED_MENUS) {
      for (const action of [
        'flip_horizontal',
        'flip_vertical',
        'add_explicit_hydrogens',
        'remove_explicit_hydrogens',
      ]) {
        await test.step(action, async () => {
          await sk.more_actions_menu('modify_all', action);
          await checkpoint(page, action);
        });
      }
      for (const format of ['sdf', 'smi', 'cxsmi', 'inchi', 'inchikey', 'pdb']) {
        await test.step(format, async () => {
          await sk.more_actions_menu('copy_all_as', format);
          await expect(await sk.clipboard_text()).not.toBe('');
        });
      }
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
      await test.step(action, async () => {
        await sk.more_actions_menu(action);
        if (reference !== null) {
          await checkpoint(page, reference);
        }
      });
    }
  });
});
