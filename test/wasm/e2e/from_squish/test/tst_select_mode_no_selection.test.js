import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker, widgetState } from '../wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
test.setTimeout(180_000);

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

async function center(sk, type, index) {
  const rect = await sk.rendered_object_rect(type, index);
  return { x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 };
}

async function drag(page, start, end, modifier = null) {
  if (modifier) await page.keyboard.down(modifier);
  try {
    await page.mouse.move(start.x, start.y, { steps: 4 });
    await page.mouse.down();
    await page.mouse.move(end.x, end.y, { steps: 12 });
    await page.mouse.up();
  } finally {
    if (modifier) await page.keyboard.up(modifier);
  }
}

test.describe('tst_select_mode_no_selection', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();

    for (const tool of ['rect_btn', 'ellipse_btn', 'lasso_btn']) {
      await sk.click_tool(tool);
      await expect((await widgetState(page, 'clear_selection_btn')).enabled).toBe(false);
      await expect((await widgetState(page, 'invert_selection_btn')).enabled).toBe(false);
    }

    for (const n of [1, 4, 7]) {
      await sk.click_atom(n, true);
      await checkpoint(page, `select_atom_${n}`);
      await sk.click_bond(n, true);
      await checkpoint(page, `select_bond_${n}`);
    }
    await sk.click_button('clear_selection');
    for (const n of [1, 4, 7]) await sk.click_atom(n, true, 'shift');
    await checkpoint(page, 'add_atoms_to_selection');
    await sk.click_button('clear_selection');
    for (const n of [1, 4, 7]) await sk.click_bond(n, true, 'shift');
    await checkpoint(page, 'add_bonds_to_selection');
    await sk.click_button('clear_selection');
    for (const n of [1, 4, 7]) { await sk.click_atom(n, true, 'shift'); await sk.click_bond(n, true, 'shift'); }
    await checkpoint(page, 'add_both_to_selection');
    for (const n of [1, 4, 7]) {
      await sk.click_bond(n, true, 'control'); await checkpoint(page, `invert_bond_${n}`);
      await sk.click_atom(n, true, 'control'); await checkpoint(page, `invert_atom_${n}`);
    }

    const oxygen = await center(sk, 'atom', 6);
    const chlorine = await center(sk, 'atom', 15);
    const nitrogen1 = await center(sk, 'atom', 1);
    const nitrogen4 = await center(sk, 'atom', 4);
    await sk.click_button('clear_selection');
    await sk.click_tool('rect_btn'); await drag(page, oxygen, chlorine); await checkpoint(page, 'mouse_drag');
    await sk.click_button('clear_selection');
    await drag(page, oxygen, chlorine); await checkpoint(page, 'shift_drag');
    await drag(page, nitrogen1, chlorine, 'Control'); await checkpoint(page, 'ctrl_drag');
    await sk.click_button('clear_selection');
    await page.mouse.dblclick(nitrogen1.x, nitrogen1.y); await checkpoint(page, 'double_click');
    await sk.click_button('clear_selection');
    await sk.click_tool('move_rotate'); await drag(page, nitrogen1, nitrogen4); await checkpoint(page, 'drag_canvas');
  });
});
