import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { hideMouseMarker } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE = 'CCC[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)CC';
test.setTimeout(120_000);

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await page.waitForTimeout(150);
  await expect(page.locator('#screen')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_edit_atom_properties', () => {
  test('atom tab basics and cancellation', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    // Preserve the live rendered-item identity across edits; modified/query
    // atom exports cannot be remapped through the original V3000 parser.
    const atomRect = await sk.rendered_object_rect('atom', sk.replay_atom_indices.get(4) || 4);
    // Keep the live screen point as the source does after its coordinate
    // rebuild: atom 4 becomes a query during this test and can no longer be
    // remapped through the original V3000 atom identity.
    const target = { type: 'screen_point', x: atomRect.x + atomRect.width / 2, y: atomRect.y + atomRect.height / 2 };
    await sk.click_atom(4, true, 'shift');

    // Source: selection context menu -> Replace with -> Allowed List.
    await sk.selection_context_menu(target, 'replace_atoms_with', 'allowed_list');
    expect((await sk.widget_state('set_as_query_rb')).checked).toBe(true);
    expect((await sk.widget_state('set_as_atom_rb')).checked).toBe(false);
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });

    await sk.selection_context_menu(target, 'edit_atom_properties');
    // Source verifies the normal Atom-tab defaults after cancelling the
    // Replace-with dialog.
    expect((await sk.widget_state('atom_element_le')).text).toBe('C');
    expect((await sk.widget_state('atom_isotope_sb')).text).toBe('');
    expect((await sk.widget_state('atom_charge_sb')).text).toBe('0');
    expect((await sk.widget_state('atom_unpaired_sb')).text).toBe('0');
    expect((await sk.widget_state('atom_stereo_combo')).text).toBe('ABS');
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });

    // Source changes every Atom-tab field, clicks Reset, and checks that the
    // original defaults are restored before cancelling the dialog.
    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({
      element: 'C', isotope: 14, charge: -2, unpaired_electrons: 2,
      enhanced_stereo: 'AND', enhanced_stereo_label: 1, click_ok: false,
    });
    await sk.edit_atom_properties({ click_ok: false, click_reset: true });
    expect((await sk.widget_state('atom_element_le')).text).toBe('C');
    expect((await sk.widget_state('atom_isotope_sb')).text).toBe('');
    expect((await sk.widget_state('atom_charge_sb')).text).toBe('0');
    expect((await sk.widget_state('atom_unpaired_sb')).text).toBe('0');
    expect((await sk.widget_state('atom_stereo_combo')).text).toBe('ABS');
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });

    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ element: 'O' });
    await checkpoint(page, 'Change_element_to_O');

    for (const [options, name] of [
      [{ unpaired_electrons: 2 }, 'Change_unpaired_electrons_to_2'],
      [{ charge: -2 }, 'Change_charge_to_minus_2'],
      [{ isotope: 14 }, 'change_isotope_to_14'],
    ]) {
      await sk.selection_context_menu(target, 'edit_atom_properties');
      await sk.edit_atom_properties(options);
      await checkpoint(page, name);
    }

    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ element: 'P', isotope: 2, charge: -5, unpaired_electrons: 4, click_ok: false });
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });
    await checkpoint(page, 'click_cancel');
  });

  test('general query types', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    await sk.click_atom(4, true);
    // Squish keeps this initial screen point after atom 4 is converted to a
    // query and uses it for each following context-menu invocation.
    const atomRect = await sk.rendered_object_rect('atom', sk.replay_atom_indices.get(4) || 4);
    const atomX = atomRect.x + atomRect.width / 2;
    const atomY = atomRect.y + atomRect.height / 2;
    await checkpoint(page, 'before_cancelling_query_tab');

    // Mirror the Squish reset sequence: alter the General query fields while
    // the dialog stays open, reset them, then accept the restored defaults.
    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({
      set_as: 'query', query_type: 'general', dropdown_type: 'Specific Element',
      element: 'F', isotope: 14, charge: 2, unpaired_electrons: 1,
      enhanced_stereo: 'AND', enhanced_stereo_label: 2,
      click_ok: false,
    });
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', click_ok: false, click_reset: true });
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', click_ok: false });
    // Source checks the General-query defaults after Reset before accepting.
    expect((await sk.widget_state('query_element_le')).text).toBe('C');
    expect((await sk.widget_state('query_isotope_sb')).text).toBe('');
    expect((await sk.widget_state('query_charge_sb')).text).toBe('');
    expect((await sk.widget_state('query_unpaired_sb')).text).toBe('');
    expect((await sk.widget_state('query_stereo_combo')).text).toBe('ABS');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general' });
    await checkpoint(page, 'after_resetting_query_tab');

    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Specific Element', element: 'F' });
    await checkpoint(page, 'change_specific_element_to_F');

    // The current dialog implementation swaps the visible input and disables
    // isotope plus unpaired electrons; charge remains editable.
    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Allowed List', click_ok: false });
    expect((await sk.widget_state('element_list_le')).enabled).toBe(true);
    expect((await sk.widget_state('element_list_le')).text).toBe('F');
    expect((await sk.widget_state('query_isotope_sb')).enabled).toBe(false);
    expect((await sk.widget_state('query_unpaired_sb')).enabled).toBe(false);
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Allowed List', element_list: 'C,N,S' });
    await checkpoint(page, 'change_allowed_list_to_CNS');

    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Not Allowed List', click_ok: false });
    expect((await sk.widget_state('element_list_le')).enabled).toBe(true);
    // Current dialog behavior restores the original specific-element value
    // when switching list polarity, rather than the prior allowed-list text.
    expect((await sk.widget_state('element_list_le')).text).toBe('F');
    expect((await sk.widget_state('query_isotope_sb')).enabled).toBe(false);
    expect((await sk.widget_state('query_unpaired_sb')).enabled).toBe(false);
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Not Allowed List', element_list: 'C,N,S' });
    await checkpoint(page, 'change_not_allowed_list_to_CNS');

    // Source accepts the Wildcard configuration (its exhaustive per-wildcard
    // loop is itself commented out), then reopens this dialog for R-Group.
    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'Wildcard', click_ok: false });
    expect((await sk.widget_state('query_isotope_sb')).enabled).toBe(true);
    // Standalone WASM keeps charge editable here; desktop MolViewer disables
    // it. This presentation difference is recorded in the parity report.
    expect((await sk.widget_state('query_charge_sb')).enabled).toBe(true);
    expect((await sk.widget_state('query_unpaired_sb')).enabled).toBe(false);
    expect((await sk.widget_state('query_stereo_combo')).enabled).toBe(false);
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general' });

    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', dropdown_type: 'R-Group', rgroup: 3 });
    await checkpoint(page, 'change_rgroup_to_3');

    // As in Squish, choosing C immediately applies it to the selected atom.
    await sk.click_tool('C');
    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({
      set_as: 'query', query_type: 'general', dropdown_type: 'Specific Element',
      element: 'F', isotope: 14, charge: -2, unpaired_electrons: 1,
    });
    await checkpoint(page, 'change_all_fields');
    await sk.selection_context_menu({ type: 'screen_point', x: atomX, y: atomY }, 'edit_atom_properties');
    await sk.edit_atom_properties({
      set_as: 'query', query_type: 'general', dropdown_type: 'Specific Element',
      element: 'C', isotope: 4, charge: -5, unpaired_electrons: 4, click_ok: false,
    });
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'general', click_ok: false, click_cancel: true });
    await checkpoint(page, 'change_all_fields_cancel');
  });

  test('advanced query number of connections', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    await sk.click_atom(4, true);
    const atomRect = await sk.rendered_object_rect('atom', sk.replay_atom_indices.get(4) || 4);
    const target = { type: 'screen_point', x: atomRect.x + atomRect.width / 2, y: atomRect.y + atomRect.height / 2 };
    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', click_ok: false });
    expect((await sk.widget_state('total_h_combo')).text).toBe('(any)');
    expect((await sk.widget_state('num_connections_sb')).text).toBe('');
    expect((await sk.widget_state('aromaticity_combo')).text).toBe('(any)');
    expect((await sk.widget_state('ring_count_combo')).text).toBe('(any)');
    expect((await sk.widget_state('ring_bond_count_combo')).text).toBe('(any)');
    expect((await sk.widget_state('smallest_ring_size_sb')).text).toBe('');
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });

    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', num_connections: 4 });
    await checkpoint(page, 'num_connections_4');
    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', total_h: '>0' });
    await checkpoint(page, 'total_h_greater_than_0');

    for (const [options, name] of [
      [{ aromaticity: 'Aliphatic (A)' }, 'aliphatic'],
      [{ aromaticity: 'Aromatic (a)' }, 'aromatic'],
      [{ aromaticity: '(any)' }, 'any_aromaticity'],
      [{ ring_count_dropdown: '>0' }, 'ring_count_greater_than_zero'],
      [{ ring_count_dropdown: 'exactly', ring_count: 4 }, 'ring_count_exactly_4'],
      [{ ring_count_dropdown: '(any)' }, 'any_ring_count'],
      [{ ring_bond_count_dropdown: 'exactly', ring_bond_count: 3 }, 'ring_bond_count_exactly_3'],
      [{ ring_bond_count_dropdown: '(any)' }, 'any_ring_bond_count'],
      [{ smallest_ring_size: 4 }, 'smallest_ring_size_4'],
    ]) {
      await sk.selection_context_menu(target, 'edit_atom_properties');
      await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', ...options });
      await checkpoint(page, name);
    }

    const allAdvancedFields = {
      total_h: '>0', num_connections: 4, aromaticity: 'Aromatic (a)',
      ring_count_dropdown: 'exactly', ring_count: 4,
      ring_bond_count_dropdown: 'exactly', ring_bond_count: 5,
      smallest_ring_size: 3,
    };
    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', ...allAdvancedFields });
    await checkpoint(page, 'set_all_fields');

    // The Squish source executes this reset/cancel path but deliberately has
    // no reference checkpoint because MolViewer does not fully reset it.
    await sk.selection_context_menu(target, 'edit_atom_properties');
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', ...allAdvancedFields, click_ok: false });
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', click_ok: false, click_reset: true });
    await sk.edit_atom_properties({ set_as: 'query', query_type: 'advanced', click_ok: false });
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });
  });

  test('atom enhanced-stereo AND and OR groups', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    const setGroups = async (type, checkpointName) => {
      await sk.import_menu('paste_in_text', SOURCE);
      await sk.map_imported_atom_indexes();
      const targets = new Map();
      for (const atom of [4, 7, 10, 13]) {
        const rect = await sk.rendered_object_rect('atom', sk.replay_atom_indices.get(atom) || atom);
        targets.set(atom, { type: 'screen_point', x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 });
      }
      for (const [atom, label] of [[4, 1], [7, 2], [10, 1], [13, 2]]) {
        await sk.click_atom(atom, true);
        await sk.selection_context_menu(targets.get(atom), 'edit_atom_properties');
        await sk.edit_atom_properties({ set_as: 'atom', enhanced_stereo: type, enhanced_stereo_label: label });
      }
      await checkpoint(page, checkpointName);
    };

    await setGroups('AND', 'test_AND_labels');
    await sk.click_button('clear');
    await setGroups('OR', 'test_OR_labels');
  });
});
