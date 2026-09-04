import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

test('selection context Copy uses a real right-click menu path', async ({ page }) => {
  const sk = new Sketcher(page);
  await sk.open();
  await sk.import_menu('paste_in_text', 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl');
  await sk.map_imported_atom_indexes();
  for (const n of [1, 2, 3, 4, 5, 6]) {
    await sk.click_bond(n, true, 'shift');
    await sk.click_atom(n, true, 'shift');
  }
  const copied = await sk.selection_context_copy_text({ type: 'atom', index: 1 });
  expect(copied).toContain('M  V30');
  const exported = await sk.selection_context_copy_as_text({ type: 'atom', index: 1 }, 'sdf');
  expect(exported).toContain('M  V30');
  await sk.selection_context_menu({ type: 'atom', index: 1 }, 'flip');
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot('flip-context.png');
  await sk.click_button('undo');
  await sk.selection_context_menu({ type: 'atom', index: 1 }, 'cut');
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot('cut-context.png');
});
