import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function checkpoint(page, name) {
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_erase_mode', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    const atoms = await sk.capture_rendered_targets('atom', [1, 3, 5, 7]);
    const bonds = await sk.capture_rendered_targets('bond', [1, 3, 5, 7]);
    await sk.click_tool('erase');
    for (const atom of [3, 5, 7]) await sk.click_rendered_target(atoms.get(atom));
    await checkpoint(page, 'delete_atoms');
    for (const bond of [3, 5, 7]) await sk.click_rendered_target(bonds.get(bond));
    await checkpoint(page, 'delete_bonds');
    await sk.drag_from_rendered_target(atoms.get(1), -600, 200);
    await checkpoint(page, 'mouse_drag_delete');

    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    const bond1 = await sk.rendered_object_rect('bond', 1);
    await sk.click_tool('triple', true);
    await sk.click_rendered_target(bond1);
    await checkpoint(page, 'change_to_triple');
    await sk.click_tool('erase');
    await sk.click_rendered_target(bond1);
    await checkpoint(page, 'erase_triple_to_double');
    await sk.click_rendered_target(bond1);
    await checkpoint(page, 'erase_double_to_single');
    await sk.click_rendered_target(bond1);
    await checkpoint(page, 'erase_single_to_none');
  });
});
