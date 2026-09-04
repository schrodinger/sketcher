import { expect, test } from '@playwright/test';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker, widgetState } from '../wrappers/sketcher_wasm.js';

const INCHI =
  'InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1';
const SMARTS = '[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]';
const SMILES =
  'CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO';
const SMALL_STRUCTURE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
// Canonical stereochemical InChI from the legacy methylphenidate.inchi
// fixture. It is a coordinate-free graph reference for every file format.
const METHYLPHENIDATE_INCHI =
  'InChI=1S/C14H19NO2/c1-17-14(16)13(11-7-3-2-4-8-11)12-9-5-6-10-15-12/h2-4,7-8,12-13,15H,5-6,9-10H2,1H3/t12-,13-/m0/s1';
const METHYLPHENIDATE_ACHIRAL_INCHI =
  'InChI=1S/C14H19NO2/c1-17-14(16)13(11-7-3-2-4-8-11)12-9-5-6-10-15-12/h2-4,7-8,12-13,15H,5-6,9-10H2,1H3';
const FILE_IMPORT_INCHI = Object.freeze({
  // PDB/ENT coordinates do not retain tetrahedral stereochemistry.
  pdb: METHYLPHENIDATE_ACHIRAL_INCHI,
  ent: METHYLPHENIDATE_ACHIRAL_INCHI,
});
const SKIP_DIALOG_STATE = process.env.PLAYWRIGHT_IMPORT_SKIP_DIALOG_STATE === '1';
const STRUCTURES_DIR = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  '../.runtime_fixtures/import_menu_structures',
);

test.setTimeout(180_000);

async function checkpoint(page, name) {
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

test.describe('tst_import_menu', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('replace-content warning', async () => {
      await sk.import_menu('paste_in_text');
      if (!SKIP_DIALOG_STATE) {
        expect(await sk.paste_in_text_status()).toBe(
          'Specified structure will <b>replace</b> Sketcher content',
        );
        expect((await widgetState(page, 'status_lbl')).styleSheet).toBe('color : #c87c00');
      }
      await sk.click_dialog_ok();
    });

    await test.step('add-content warning', async () => {
      await sk.import_menu('replace_current_content');
      await sk.import_menu('paste_in_text');
      if (!SKIP_DIALOG_STATE) {
        expect(await sk.paste_in_text_status()).toBe(
          'Specified structure will be added to Sketcher content',
        );
        expect((await widgetState(page, 'status_lbl')).styleSheet).toBe('color : gray');
      }
      await sk.click_dialog_ok();
      await sk.import_menu('replace_current_content');
    });

    for (const [name, text, format] of [
      ['import_inchi_text', INCHI, 'Autodetect (InChI)'],
      ['import_smarts_text', SMARTS, 'Autodetect (SMARTS)'],
      ['import_smiles_text', SMILES, 'Autodetect (SMILES)'],
    ]) {
      await test.step(name, async () => {
        await sk.import_menu('paste_in_text', text, false);
        expect(await sk.paste_in_text_format()).toBe(format);
        await sk.click_dialog_ok();
        await checkpoint(page, name);
        if (name === 'import_inchi_text') {
          // Exact InChI comparison includes the source molecule's tetrahedral
          // stereochemistry, beyond the visual wedge/dash checkpoint.
          await expect(await sk.copy_all_as_text('inchi')).toBe(INCHI);
        }
      });
    }

    await test.step('import_sdf_text', async () => {
      const sdfText = await sk.copy_all_as_text('sdf');
      expect(sdfText).toContain('V3000');
      await sk.import_menu('paste_in_text', sdfText, false);
      expect(await sk.paste_in_text_format()).toBe('Autodetect (MDL SD)');
      await sk.click_dialog_ok();
      await checkpoint(page, 'import_sdf_text');
    });

    await test.step('replace and add small structures', async () => {
      await sk.import_menu('paste_in_text', SMALL_STRUCTURE);
      await checkpoint(page, 'import_small_structure');
      await sk.import_menu('replace_current_content');
      for (const name of [
        'import_small_structure_no_replace_1',
        'import_small_structure_no_replace_2',
        'import_small_structure_no_replace_3',
      ]) {
        await sk.import_menu('paste_in_text', SMALL_STRUCTURE);
        await checkpoint(page, name);
      }
      await sk.click_button('cleanup');
      await checkpoint(page, 'import_small_structure_no_replace_cleanup');
      await sk.import_menu('replace_current_content');
    });

    for (const extension of [
      'pdb',
      'smi',
      'smiles',
      'cxsmiles',
      'sdf',
      'mdl',
      'ent',
      'inchi',
      'mol',
      'sd',
      'cxsmi',
    ]) {
      await test.step(`import_from_file_${extension}`, async () => {
        await sk.import_menu(
          'import_from_file',
          path.join(STRUCTURES_DIR, `methylphenidate.${extension}`),
        );
        await checkpoint(page, `import_from_file_${extension}`);
        await expect(await sk.copy_all_as_text('inchi')).toBe(
          FILE_IMPORT_INCHI[extension] || METHYLPHENIDATE_INCHI,
        );
      });
    }
  });
});
