import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  clearSketcher,
  tryImport,
  getExportedSmiles,
  isSketcherEmpty,
} from './e2e_helpers.js';

// Ported from `suite_2D_sketcher_new/tst_import_menu`. The original drives the
// sketcher's Import dialog (paste-in-text box, file picker, "replace content"
// checkbox) inside Maestro and compares pixel screenshots. The standalone WASM
// sketcher has no import dialog — structures come in through
// Module.sketcher_import_text — so we preserve the testable intent: every
// supported text/file format imports to the correct structure, and the
// 200-atom guard rejects oversized input.

const FIXTURES_DIR = join(
  dirname(fileURLToPath(import.meta.url)),
  'fixtures',
  'import_menu_structures',
);

// Methylphenidate, supplied by the squish suite in 11 formats. The lossless
// formats must all canonicalize to one SMILES; PDB/ENT (added hydrogens,
// coordinate-derived bonds) and InChI (stereo perception) may legitimately
// differ, so they are only required to import non-empty.
const LOSSLESS = ['smi', 'smiles', 'cxsmi', 'cxsmiles', 'sdf', 'sd', 'mol', 'mdl'];
const LOSSY = ['pdb', 'ent', 'inchi'];

const fixture = (ext) => readFileSync(join(FIXTURES_DIR, `methylphenidate.${ext}`), 'utf8');

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Import formats', () => {
  test('lossless formats all canonicalize to the same SMILES', async ({ page }) => {
    // Establish the reference from the plain SMILES fixture, then require every
    // other lossless format to match it (rather than hard-coding the sketcher's
    // canonical form).
    await clearSketcher(page);
    await tryImport(page, fixture('smi'));
    const reference = await getExportedSmiles(page);
    expect(reference, 'reference SMILES should be non-empty').toBeTruthy();

    for (const ext of LOSSLESS) {
      await clearSketcher(page);
      const result = await tryImport(page, fixture(ext));
      expect(result.ok, `.${ext} import threw: ${result.error}`).toBe(true);
      expect(
        await getExportedSmiles(page),
        `.${ext} canonical SMILES should match the .smi reference`,
      ).toBe(reference);
    }
  });

  for (const ext of LOSSY) {
    test(`lossy format .${ext} imports non-empty`, async ({ page }) => {
      await clearSketcher(page);
      const result = await tryImport(page, fixture(ext));
      expect(result.ok, `.${ext} import threw: ${result.error}`).toBe(true);
      expect(result.empty, `.${ext} produced an empty canvas`).toBe(false);
    });
  }

  // The InChI / SMARTS / SMILES sample strings exercised by the original test.
  test('InChI string imports', async ({ page }) => {
    await clearSketcher(page);
    const inchi =
      'InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1';
    const result = await tryImport(page, inchi);
    expect(result.ok, `InChI import threw: ${result.error}`).toBe(true);
    expect(result.empty).toBe(false);
  });

  test('SMARTS query string imports', async ({ page }) => {
    await clearSketcher(page);
    const smarts = '[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]';
    const result = await tryImport(page, smarts);
    expect(result.ok, `SMARTS import threw: ${result.error}`).toBe(true);
    expect(result.empty).toBe(false);
  });

  test('complex SMILES (cephalostatin-1) imports', async ({ page }) => {
    await clearSketcher(page);
    const smiles =
      'CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO';
    const result = await tryImport(page, smiles);
    expect(result.ok, `SMILES import threw: ${result.error}`).toBe(true);
    expect(result.empty).toBe(false);
  });

  test('oversized structure is rejected by the 200-atom limit', async ({ page }) => {
    // A >200-atom polysaccharide SMILES; the original asserted the sketcher
    // surfaced "too many atoms present" and imported nothing.
    const longText =
      'OC[C@@]1([H])O[C@@](OC[C@@]2([H])O[C@@](OC[C@@]3([H])O[C@@](OC[C@@]4([H])O[C@@](OC[C@@]5([H])O[C@@](OC[C@@]6([H])O[C@@](OC[C@@]7([H])O[C@@](OC[C@@]8([H])O[C@@](OC[C@@]9([H])O[C@@](OC[C@@]%10([H])O[C@]([H])(OC[C@@]%11([H])O[C@@](OC[C@@]%12([H])O[C@@](OC[C@@]%13([H])O[C@@](OC[C@@]%14([H])O[C@@](OC[C@@]%15([H])O[C@@](OC[C@@]%16([H])O[C@@](OC[C@@]%17([H])O[C@@](OC[C@@]%18([H])O[C@@](OC[C@@]%19([H])O[C@@](OC[C@@]%20([H])O[C@@](O)([H])[C@]([H])(O)[C@](O)([H])[C@]%20([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%19([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%18([H])O)([H])[C@]([H])(O)[C@@]([H])(O[C@H]%21O[C@]([H])(CO)[C@@]([H])(O)[C@@](O)([H])[C@]%21(O)[H])[C@]%17([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%16([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%15([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%14([H])O)([H])[C@]([H])(O[C@H]%22O[C@]([H])(CO)[C@@]([H])(O)[C@@](O)([H])[C@]%22(O)[H])[C@](O)([H])[C@]%13([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%12([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]%11([H])O)[C@]([H])(O)[C@](O)([H])[C@]%10([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]9([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]8([H])O)([H])[C@]([H])(O)[C@@]([H])(O[C@H]%23O[C@]([H])(CO)[C@@]([H])(O)[C@@](O)([H])[C@]%23(O)[H])[C@]7([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]6([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]5([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]4([H])O)([H])[C@]([H])(O[C@H]%24O[C@]([H])(CO)[C@@]([H])(O)[C@@](O)([H])[C@]%24(O)[H])[C@](O)([H])[C@]3([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]2([H])O)([H])[C@]([H])(O)[C@](O)([H])[C@]1([H])O';

    await clearSketcher(page);
    await tryImport(page, longText);
    // Whether the limit surfaces as a thrown error or a silently-ignored
    // import, the canvas must remain empty.
    expect(await isSketcherEmpty(page), 'oversized structure should not import').toBe(true);
  });
});
