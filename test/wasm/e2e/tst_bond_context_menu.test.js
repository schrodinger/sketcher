import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { hideMouseMarker } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_bond_context_menu', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source: import, retrieve coordinates, clear, then rebuild. The browser
    // wrapper records the imported atom-to-rendered-item mapping instead.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    const target = { type: 'bond', index: 4 };

    await test.step('flip_substituent', async () => {
      await sk.click_bond(4, true);
      await sk.selection_context_menu(target, 'flip_substituent');
      await checkpoint(page, 'flip_substituent');
    });

    // Source: select bonds one through six before testing the remaining
    // context-menu actions.
    for (const index of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(index, true, 'shift');
    }

    for (const action of ['double', 'single', 'aromatic', 'up', 'down']) {
      await test.step(action, async () => {
        await sk.selection_context_menu(target, action);
        await checkpoint(page, action);
      });
    }

    await test.step('delete', async () => {
      // Change back to single bonds before deletion, exactly as in Squish.
      await sk.selection_context_menu(target, 'single');
      await sk.selection_context_menu(target, 'delete');
      await checkpoint(page, 'delete');
      await sk.click_button('undo');
    });

    for (const action of ['Coordinate', 'Zero Order', 'Single Up/Down', 'Double Cis/Trans']) {
      await test.step(`other_type: ${action}`, async () => {
        await sk.selection_context_menu(target, 'other_type', action);
        await checkpoint(page, action.replace('/', '_'));
      });
    }

    for (const action of ['Any', 'Single/Double', 'Double/Aromatic', 'Single/Aromatic']) {
      await test.step(`query: ${action}`, async () => {
        await sk.selection_context_menu(target, 'query', action);
        await checkpoint(page, action.replace('/', '_'));
      });
    }

    // Source: return to single bonds before topology checks.
    await sk.selection_context_menu(target, 'single');

    for (const action of ['In Ring', 'Not In a Ring', 'Either']) {
      await test.step(`topology: ${action}`, async () => {
        await sk.selection_context_menu(target, 'topology', action);
        await checkpoint(page, action);
      });
    }
  });
});
