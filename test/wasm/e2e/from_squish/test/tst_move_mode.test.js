import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_move_mode', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source: import, obtain source coordinates, clear, then rebuild.  Keep
    // the imported molecule here—the geometry bridge gives the same stable
    // targets without adding an otherwise unrelated construction sequence.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();

    // Select atoms and bonds 1, 4, and 7 before selecting Move/Rotate.
    for (const n of [1, 4, 7]) {
      await sk.click_bond(n, true, 'shift');
      await sk.click_atom(n, true, 'shift');
    }
    await sk.click_tool('move_rotate');
    await checkpoint(page, 'maintain_previous_selection');

    // Drag selection from source atom 4 by (+50, +50) in Squish coordinates.
    const atom4 = sk.replay_atom_coordinates.get(4);
    const atom3 = sk.replay_atom_coordinates.get(3);
    await sk.mouse_drag(atom4.x, atom4.y, 50, 50);
    await checkpoint(page, 'drag_selection');
    await sk.click_button('undo');

    // Rotation-circle gesture: source atom 4 + (130, 0) to source atom 3.
    await sk.mouse_drag(atom4.x + 130, atom4.y, atom3.x - (atom4.x + 130), atom3.y - atom4.y);
    await checkpoint(page, 'rotate_selection');
    await sk.click_button('undo');

    // Preserve Squish's explicit background drag workaround.
    await sk.drag_canvas_point(10, 10, 110, 110);
    await checkpoint(page, 'drag_background');

    // Uncheck Replace Current Content and import a second structure, as in
    // the source test.  Live item 18 is the lower nitrogen selected there.
    await sk.click_button('clear_selection');
    await sk.click_tool('rect_btn');
    await sk.import_menu('replace_current_content');
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.click_button('cleanup');
    // Squish obtains a new structure dictionary here before replaying the
    // second structure.  Refresh the live coordinate map instead of doing
    // its clear/rebuild workaround.
    await sk.map_imported_atom_indexes();
    const atom18 = await sk.rendered_object_rect('atom', 18);
    await sk.click_tool('move_rotate');

    // The source drags from atom 18 + (30, -40) to source position (50, -50).
    // Use the current rendered atom to retain that human mouse gesture despite
    // standalone canvas scaling differing from desktop MolViewer.
    await sk.drag_rendered_target_to_source_point(atom18, 30, -40, 50, -50);
    await checkpoint(page, 'drag_rotate_tool_multiple_structures');
    await sk.click_button('undo');
    await sk.click_tool('rect_btn');
    await sk.click_tool('move_rotate');
    await sk.drag_canvas_point(10, 10, 110, 110);
    await checkpoint(page, 'drag_background_multiple_structures');
  });
});
