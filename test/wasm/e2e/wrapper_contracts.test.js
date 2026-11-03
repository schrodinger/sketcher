import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { setClipboardText } from './sketcher/wrappers/sketcher_wasm.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
test.setTimeout(90_000);

function v3000Counts(text) {
  const match = text.match(/^M  V30 COUNTS (\d+) (\d+) /m);
  return match ? { atoms: Number(match[1]), bonds: Number(match[2]) } : null;
}

async function visibleQtPopupCount(page) {
  return page.locator('[id^="qt-window-"]').evaluateAll(
    (elements) =>
      elements.slice(1).filter((element) => {
        const style = getComputedStyle(element);
        const bounds = element.getBoundingClientRect();
        return (
          style.display !== 'none' &&
          style.visibility !== 'hidden' &&
          bounds.width > 0 &&
          bounds.height > 0
        );
      }).length,
  );
}

test.describe('Sketcher wrapper contracts', () => {
  test('Import menu -> Paste in Text accepts text through the visible dialog', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await expect(await sk.copy_all_as_text('smi')).toBe(SOURCE_STRUCTURE);
  });

  test('Select tool dropdown switches between lasso, ellipse, and rectangle', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);

    for (const tool of ['lasso_btn', 'ellipse_btn', 'rect_btn']) {
      await sk.click_tool(tool);
      expect((await sk.widget_state(tool)).checked).toBe(true);
      await expect.poll(() => visibleQtPopupCount(page)).toBe(0);
    }
  });

  test('selection after Clear and keyboard shortcuts preserves the first bond click', async ({
    page,
  }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', 'CCO');
    await sk.click_tool('rect_btn');
    await sk.click_button('clear');
    await sk.import_menu('paste_in_text', 'CCO');
    await sk.type_text('sketcher_area', '<Ctrl+A>');
    await sk.type_text('sketcher_area', '<Ctrl+D>');
    expect((await sk.widget_state('select_tool_btn')).checked).toBe(true);

    await sk.click_tool('rect_btn');
    await expect.poll(() => visibleQtPopupCount(page)).toBe(0);
    await sk.click_bond(1, true, 'shift');
    const sentinel = '__playwright_waiting_for_bond_copy__';
    await setClipboardText(page, sentinel);
    await sk.type_text('sketcher_area', '<Ctrl+C>');
    // Qt's copy event can publish the cached import text before Sketcher's
    // asynchronous clipboard write finishes. Wait for the expected payload,
    // since merely replacing the sentinel can still return that stale text.
    await expect
      .poll(async () => v3000Counts(await sk.clipboard_text()))
      .toEqual({
        atoms: 2,
        bonds: 1,
      });
  });

  test('combo choices leave atom properties open and can be cancelled', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', 'C[C@H](F)Cl');
    const original = await sk.copy_all_as_text('smi');
    await sk.map_imported_atom_indexes();
    await sk.click_atom(2, true);

    for (const enhanced_stereo of ['AND', 'OR']) {
      await sk.selection_context_menu({ type: 'atom', index: 2 }, 'edit_atom_properties');
      await sk.edit_atom_properties({ enhanced_stereo, click_ok: false });
      expect((await sk.widget_state('atom_stereo_combo')).text).toBe(enhanced_stereo);
      await sk.edit_atom_properties({ click_ok: false, click_cancel: true });
    }

    await sk.selection_context_menu({ type: 'atom', index: 2 }, 'edit_atom_properties');
    await sk.edit_atom_properties({
      set_as: 'query',
      query_type: 'advanced',
      total_h: '>0',
      aromaticity: 'Aromatic (a)',
      ring_count_dropdown: 'exactly',
      ring_count: 3,
      click_ok: false,
    });
    expect((await sk.widget_state('total_h_combo')).text).toBe('>0');
    expect((await sk.widget_state('aromaticity_combo')).text).toBe('Aromatic (a)');
    await sk.edit_atom_properties({ click_ok: false, click_cancel: true });
    await sk.click_button('clear_selection');
    expect(await sk.copy_all_as_text('smi')).toBe(original);
  });

  test('reset_state restores the Squish sticky-tool defaults', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await sk.click_tool('ellipse_btn');
    await sk.click_tool('XH');
    await sk.click_tool('triple');

    await sk.reset_state();

    expect(sk.current_buttons).toMatchObject({
      tool: 'C',
      wildcard: 'A',
      periodic_table: 'Si',
      stereo: 'down',
      bond_order: 'double',
      bond_query: 'aromatic',
      'r-group': 'r+',
      reaction: 'rxn_arrow',
      select: 'rect_btn',
    });
    expect((await sk.widget_state('c_btn')).checked).toBe(true);
  });

  test('Space preserves Select mode for a following lasso selection', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await sk.map_imported_atom_indexes();

    // Match the source shortcut test's precondition: Space is exercised after
    // Backspace has acted on an active selection, not immediately after import.
    await sk.click_atom(1, true, 'shift');
    await sk.click_bond(1, true, 'shift');
    await sk.type_text('sketcher_area', '<Backspace>');

    await sk.click_tool('C');
    await sk.type_text('sketcher_area', ' ');
    expect((await sk.widget_state('select_tool_btn')).checked).toBe(true);

    await sk.click_tool('lasso_btn');
    await sk.click_tool('C');
    await sk.type_text('sketcher_area', ' ');
    expect((await sk.widget_state('select_tool_btn')).checked).toBe(true);
  });

  test('Bond context Topology accepts the source Not In a Ring label', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await sk.map_imported_atom_indexes();
    // The source makes this a multi-bond selection before exposing its
    // Modify Bonds submenu, so preserve that UI precondition here.
    for (const index of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(index, true, 'shift');
    }

    await sk.selection_context_menu({ type: 'bond', index: 4 }, 'topology', 'Not In a Ring');
  });

  test('deep context-menu action dismisses its Qt/WASM popup canvases', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await sk.map_imported_atom_indexes();
    // Preserve the source precondition for the nested Modify Atoms menu.
    for (const index of [1, 2, 3, 4, 5, 6]) {
      await sk.click_bond(index, true, 'shift');
      await sk.click_atom(index, true, 'shift');
    }

    await sk.selection_context_menu(
      { type: 'atom', index: 1 },
      'modify_atoms',
      'replace_atoms_with',
      'wildcard',
      'A',
    );
    await expect.poll(() => visibleQtPopupCount(page)).toBe(0);
  });

  test('Bracket Subgroup dialog creates the source SRU/HT S-group', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await sk.map_imported_atom_indexes();

    await sk.add_sgroup('SRU', [4, 5, 6, 7], [4, 5, 6], 'HT', '');
    const sdf = await sk.copy_all_as_text('sdf');
    expect(sdf).toContain('M  V30 BEGIN SGROUP');
    expect(sdf).toContain(' SRU ');
    expect(sdf).toContain('CONNECT=HT');
  });

  test('Help menu opens and closes both standalone dialogs', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await sk.help_menu('getting_started');
    await sk.click_button('ok_btn');
    await sk.help_menu('about_sketcher');
    await sk.click_button('about_sketcher_ok');
  });

  test('More Actions submenu representative paths', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);

    await test.step('Copy All As -> SDF', async () => {
      await expect(v3000Counts(await sk.copy_all_as_text('sdf'))).toEqual({ atoms: 15, bonds: 15 });
    });

    await test.step('Modify All -> Add and Remove Explicit Hydrogens', async () => {
      await sk.more_actions_menu('modify_all', 'add_explicit_hydrogens');
      const withHydrogens = v3000Counts(await sk.copy_all_as_text('sdf'));
      await expect(withHydrogens?.atoms).toBeGreaterThan(15);
      await expect(withHydrogens?.bonds).toBeGreaterThan(15);

      await sk.more_actions_menu('modify_all', 'remove_explicit_hydrogens');
      await expect(v3000Counts(await sk.copy_all_as_text('sdf'))).toEqual({ atoms: 15, bonds: 15 });
    });
  });

  test('Export menu selects SDF through the visible File Export dialog', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);

    const download = await sk.export_menu({ filename: 'contract_export', format: 'SDF' });
    expect(download.filename).toBe('contract_export.sdf');
    expect(Buffer.from(download.contentBase64, 'base64').toString('utf8')).toContain('V3000');
  });

  test('Save Image exposes the source size and transparency controls', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);

    const download = await sk.save_image({
      filename: 'contract_image',
      format: 'PNG',
      width: 420,
      height: 690,
      transparent: true,
    });
    expect(download.filename).toBe('contract_image.png');
    expect(Buffer.from(download.contentBase64, 'base64').subarray(0, 8)).toEqual(
      Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    );
  });
});
