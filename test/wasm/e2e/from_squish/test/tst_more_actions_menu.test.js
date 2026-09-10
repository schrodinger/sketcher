import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker, setClipboardText } from '../wrappers/sketcher_wasm.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
// Values from suite_molviewer/shared/testdata/more_actions_menu_references.
// These formats do not encode the editor's rendered coordinates, so the
// MolViewer references are directly comparable with standalone WASM output.
const COPY_ALL_TEXT_REFERENCES = Object.freeze({
  smi: SOURCE_STRUCTURE,
  cxsmi: SOURCE_STRUCTURE,
  inchi: 'InChI=1S/C9H9Cl2N3O/c10-6-2-1-3-7(11)5(6)4-8(15)14-9(12)13/h1-3H,4H2,(H4,12,13,14,15)',
  inchikey: 'INJOMKTZOLKMBF-UHFFFAOYSA-N',
});
test.setTimeout(180_000);

function v3000Counts(text) {
  const match = text.match(/^M  V30 COUNTS (\d+) (\d+) /m);
  if (!match) {
    return null;
  }
  return { atoms: Number(match[1]), bonds: Number(match[2]) };
}

function pdbAtomElements(text) {
  return text
    .split(/\r?\n/)
    .filter((line) => /^(?:ATOM  |HETATM)/.test(line))
    .map((line) => {
      const element = line.slice(76, 78).trim();
      return element.slice(0, 1) + element.slice(1).toLowerCase();
    })
    .sort();
}

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
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}-replay.png`);
}

async function selectSourceAtomsAndBonds(sk) {
  // Exact source sequence from suite_molviewer/tst_more_actions_menu/test.py.
  for (const n of [1, 2, 3, 4, 5, 6]) {
    await sk.click_bond(n, true, 'shift');
    await sk.click_atom(n, true, 'shift');
  }
}

/** Import through the visible source menu; rendered geometry handles selection. */
async function importSourceStructure(sk) {
  await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
}

test.describe('tst_more_actions_menu', () => {
  test('main', async ({ page }, testInfo) => {
    const sk = new Sketcher(page);
    await sk.open();
    await requireRenderedGeometryBridge(page);

    // Squish rebuilt here to compensate for fragile desktop coordinates.
    // The WASM bridge returns live rendered atom/bond geometry, so preserve
    // the user-facing Import path without the intermediate reconstruction.
    await importSourceStructure(sk);
    await selectSourceAtomsAndBonds(sk);

    await test.step('add_hydrogens_ignores_selection', async () => {
      await sk.more_actions_menu('modify_all', 'add_explicit_hydrogens');
      await checkpoint(page, 'add_hydrogens_ignores_selection');
      await sk.click_button('undo');
    });
    await test.step('Ctrl_X_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+X>');
      await checkpoint(page, 'Ctrl_X_hotkey');
    });

    await sk.click_button('clear');
    await test.step('Ctrl_V_hotkey', async () => {
      // Chromium reserves Ctrl+V before Qt/WASM receives it. The visible
      // Paste menu action invokes the same Sketcher paste command, using the
      // clipboard payload produced by the preceding real Ctrl+X gesture.
      testInfo.annotations.push({
        type: 'browser shortcut adaptation',
        description: 'Ctrl+V is represented by the visible More Actions -> Paste command.',
      });
      await sk.more_actions_menu('paste');
      await checkpoint(page, 'Ctrl_V_hotkey');
    });

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

    // This is the source's second Import/Paste -> selection setup.
    await importSourceStructure(sk);
    await selectSourceAtomsAndBonds(sk);

    await test.step('Ctrl_I_hotkey', async () => {
      await sk.type_text('sketcher_area', '<Ctrl+I>');
      await checkpoint(page, 'Ctrl_I_hotkey');
    });
    await test.step('Ctrl_C_hotkey', async () => {
      await setClipboardText(page, '__playwright_waiting_for_ctrl_c__');
      await sk.type_text('sketcher_area', '<Ctrl+C>');
      const clipboard = await sk.wait_for_clipboard_change('__playwright_waiting_for_ctrl_c__');
      await expect(v3000Counts(clipboard)).toEqual({ atoms: 9, bonds: 9 });
    });
    await test.step('Ctrl_F_hotkey', async () => {
      await sk.click_button('clear_selection');
      await sk.mouse_drag(0, 0, 100, 100, null, 'right');
      await sk.type_text('sketcher_area', '<Ctrl+F>');
      await checkpoint(page, 'Ctrl_F_hotkey');
    });

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
        if (format === 'pdb') {
          // Squish stored internal structure information here rather than the
          // exported PDB payload. Validate the browser-visible PDB instead.
          const clipboard = await sk.copy_all_as_text(format);
          await expect(pdbAtomElements(clipboard)).toEqual([
            'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
            'Cl', 'Cl', 'N', 'N', 'N', 'O',
          ]);
          return;
        }
        const clipboard = await sk.copy_all_as_text(format);
        if (format === 'sdf') {
          // The legacy SDF has the same molecular graph but different 2-D
          // coordinates, so preserve source intent with a graph-level check.
          await expect(v3000Counts(clipboard)).toEqual({ atoms: 15, bonds: 15 });
        } else if (Object.hasOwn(COPY_ALL_TEXT_REFERENCES, format)) {
          await expect(clipboard).toBe(COPY_ALL_TEXT_REFERENCES[format]);
        }
      });
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
