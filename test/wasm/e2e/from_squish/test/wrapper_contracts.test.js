import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';

const SOURCE_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
test.setTimeout(90_000);

function v3000Counts(text) {
  const match = text.match(/^M  V30 COUNTS (\d+) (\d+) /m);
  return match ? { atoms: Number(match[1]), bonds: Number(match[2]) } : null;
}

test.describe('Sketcher wrapper contracts', () => {
  test('Import menu -> Paste in Text accepts text through the visible dialog', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);
    await expect(await sk.copy_all_as_text('smi')).toBe(SOURCE_STRUCTURE);
  });

  test('Select tool dropdown exposes lasso and ellipse choices', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();
    await sk.import_menu('paste_in_text', SOURCE_STRUCTURE);

    await sk.click_tool('lasso_btn');
    await sk.click_tool('ellipse_btn');
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

    await sk.selection_context_menu(
      { type: 'bond', index: 4 },
      'topology',
      'Not In a Ring',
    );
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
});
