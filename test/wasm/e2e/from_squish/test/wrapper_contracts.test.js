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
