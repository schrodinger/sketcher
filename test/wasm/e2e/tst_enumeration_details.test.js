import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { hideMouseMarker } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
test.setTimeout(180_000);

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

async function importAndMap(sk) {
  await sk.import_menu('paste_in_text', SOURCE);
  // The desktop source obtains coordinates, clears, and redraws the structure.
  // Retain the imported molecule and map those source atom numbers to its live
  // browser geometry; this preserves the source actions without fragile pixels.
  await sk.map_imported_atom_indexes();
}

test.describe('tst_enumeration_details', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await importAndMap(sk);

    await sk.click_tool('r+');
    for (const [atom, name] of [[1, 'add_r1'], [2, 'add_r2'], [3, 'add_r3']]) {
      await sk.click_atom(atom);
      await checkpoint(page, name);
    }

    await sk.click_tool('r1');
    await sk.click_atom(1);
    await checkpoint(page, 'add_r1_to_r1');
    await sk.click_atom(1);
    await checkpoint(page, 'add_r1_to_r1_again');

    await sk.click_tool('r9');
    await sk.click_atom(4);
    await checkpoint(page, 'add_r9');

    await sk.click_tool('r+');
    for (const atom of [5, 6, 7, 8, 9, 10, 11, 12]) await sk.click_atom(atom);
    await checkpoint(page, 'fill_in_r_groups');

    await sk.click_tool('r5');
    for (const atom of [8, 9, 10]) await sk.click_atom(atom);
    await checkpoint(page, 'change_to_r5');

    await sk.reset_state();
    await importAndMap(sk);

    await sk.click_tool('attachment_point');
    await sk.click_atom(5);
    await checkpoint(page, 'add_attachment_point_draw_mode');
    await sk.click_atom(5);
    await checkpoint(page, 'add_attachment_point_draw_mode_again');

    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', 'CCC.CCCCC>>CCCCCCC');
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');

    await sk.click_button('clear');
    await sk.click_tool('atom_chain');
    await sk.mouse_drag(-300, 0, 100, 0);
    await sk.click_tool('rxn_plus');
    await sk.click_sketcher(-175, 0);
    await sk.click_tool('C');
    await sk.mouse_drag(-150, 0, 75, 0);
    await sk.click_tool('rxn_arrow');
    await sk.click_sketcher(0, 0);
    await sk.click_tool('atom_chain');
    await sk.mouse_drag(50, 0, 175, 0);
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');

    await sk.click_tool('add_mapping');
    await sk.mouse_drag(-300, 0, 350, 0);
    await sk.mouse_drag(-300, 0, 435, 0);
    await sk.mouse_drag(-150, 0, 370, 0);
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');

    await sk.click_tool('remove_mapping');
    await sk.click_sketcher(50, 0);
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');
    await sk.click_sketcher(220, 0);
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');
    await sk.click_sketcher(135, 0);
    await expect(await sk.copy_all_as_text('smi')).toContain('>>');
  });
});
