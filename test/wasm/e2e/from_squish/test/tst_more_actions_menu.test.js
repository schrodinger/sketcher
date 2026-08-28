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
});
