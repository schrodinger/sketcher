import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
const PREFERENCES = {
  atom_font_size: 16,
  bond_line_width: 2.5,
  label_carbons: true,
  terminal_only: false,
  all_carbons: true,
  with_abs: true,
  include_undefined_centers: false,
  color_mode: 'dark',
};

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_configure_view_menu', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source: import, retrieve atom coordinates, clear, and rebuild. The
    // browser wrapper maps the imported atom indices before each interaction.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();

    await test.step('change_bond_to_up', async () => {
      await sk.click_tool('up');
      await sk.click_bond(5);
      await checkpoint(page, 'change_bond_to_up');
    });

    await test.step('stereo_labels_off', async () => {
      await sk.configure_view_menu('stereo_labels');
      await checkpoint(page, 'stereo_labels_off');
    });

    await test.step('stereo_labels_on', async () => {
      await sk.configure_view_menu('stereo_labels');
      await checkpoint(page, 'stereo_labels_on');
    });

    await test.step('heteroatoms_with_color', async () => {
      for (const [index, element] of ['N', 'O', 'P', 'F', 'Cl', 'Br', 'I'].entries()) {
        await sk.click_tool(element);
        await sk.click_atom(index + 1);
      }
      await checkpoint(page, 'heteroatoms_with_color');
    });

    await test.step('heteroatom_colors_off', async () => {
      await sk.configure_view_menu('heteroatom_colors');
      await checkpoint(page, 'heteroatom_colors_off');
    });

    await test.step('heteroatom_colors_on', async () => {
      await sk.configure_view_menu('heteroatom_colors');
      await checkpoint(page, 'heteroatom_colors_on');
    });

    await test.step('valence_errors_off', async () => {
      await sk.configure_view_menu('valence_errors');
      await checkpoint(page, 'valence_errors_off');
    });

    await test.step('valence_errors_on', async () => {
      await sk.configure_view_menu('valence_errors');
      await checkpoint(page, 'valence_errors_on');
    });

    // The source changes the chlorine back before opening preferences.
    await sk.click_tool('C');
    await sk.click_atom(5);

    await test.step('check_default_preferences', async () => {
      await sk.configure_view_menu('preferences', { close: false });
      expect((await sk.widget_state('m_atom_font_size_sb')).text).toBe('18');
      expect((await sk.widget_state('m_bond_line_width_sb')).text).toBe('2.40');
      expect((await sk.widget_state('m_label_carbons_cb')).checked).toBe(false);
      expect((await sk.widget_state('m_show_stereo_cb')).checked).toBe(true);
      expect((await sk.widget_state('m_color_heteroatoms_cb')).checked).toBe(true);
      await sk.click_button('m_close_btn');
    });

    await test.step('change_preferences', async () => {
      await sk.configure_view_menu('preferences', PREFERENCES);
      await checkpoint(page, 'change_preferences');
    });

    await test.step('reset_preferences', async () => {
      await sk.configure_view_menu('preferences', { reset: true });
      await checkpoint(page, 'reset_preferences');
    });

    await test.step('close_preferences', async () => {
      await sk.configure_view_menu('preferences', {
        ...PREFERENCES,
        atom_font_size: 8,
        bond_line_width: 4.5,
        close: false,
      });
      expect((await sk.widget_state('m_atom_font_size_sb')).text).toBe('8');
      expect((await sk.widget_state('m_bond_line_width_sb')).text).toBe('4.50');
      await sk.click_button('m_close_btn');
      await checkpoint(page, 'close_preferences');
    });
  });
});
