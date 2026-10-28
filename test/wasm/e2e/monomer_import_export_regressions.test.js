import { expect, test } from '@playwright/test';
import path from 'node:path';
import {
  FIXTURES,
  clear,
  clickWidget,
  closeDialog,
  expectDialog,
  getExportedHelm,
  getExportedSmiles,
  openSketcher,
  pasteText,
  scene,
  screenshotView,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';
import { loadDatabase, openMonomerSketcher } from './monomer_helpers.js';
import { clickPopupRow } from './sketcher/wrappers/sketcher_wasm.js';

test.setTimeout(120_000);

test('SKETCH-1848: the original 3D SDF displays R and retains stereochemistry through cleanup and SMILES', async ({
  page,
}) => {
  const sk = await openSketcher(page);
  await sk.import_menu('import_from_file', path.join(FIXTURES, 'mol2-1.sdf'));
  await waitAtoms(page, 26);
  const chiralLabels = async () =>
    (await scene(page)).atoms
      .filter((atom) => atom.chiralityText)
      .map((atom) => atom.chiralityText);
  // Read AtomItem's actual painted label, not just the molecule's CIP property.
  await expect.poll(chiralLabels).toEqual(['(R)']);
  await clickWidget(page, 'fit_btn');
  await page.mouse.move(1100, 40);
  await screenshotView(page, '3d-import-R-label.png');
  const smiles = await sk.copy_all_as_text('smi');
  expect(smiles).toContain('@');
  await clickWidget(page, 'cleanup_btn');
  expect(await getExportedSmiles(page)).toBe(smiles);
  expect(await chiralLabels()).toEqual(['(R)']);
  await clear(page);
  await pasteText(page, smiles);
  expect(await chiralLabels()).toEqual(['(R)']);
  expect(await getExportedSmiles(page)).toBe(smiles);
});

test.describe('Monomer import validation — SKETCH-2862', () => {
  for (const [name, input, error] of [
    [
      'missing peptide',
      'PEPTIDE1{[THIS_IS_NOT_A_REAL_RESIDUE_NAME]}$$$$V2.0',
      /Peptide monomer THIS_IS_NOT_A_REAL_RESIDUE_NAME not found in monomer database/,
    ],
    [
      'missing base',
      'RNA1{R([missingMonomer])P}$$$$V2.0',
      /Nucleic acid monomer missingMonomer not found in monomer database/,
    ],
    [
      'missing CHEM',
      'CHEM1{[missingMonomer]}$$$$V2.0',
      /CHEM monomer missingMonomer not found in monomer database/,
    ],
    [
      'invalid inline SMILES',
      'PEPTIDE1{[C1CC]}$$$$V2.0',
      /SMILES|parse|unrecognized|format|not found/i,
    ],
  ]) {
    test(`${name} shows an error and leaves existing content intact in Add mode`, async ({
      page,
    }) => {
      await openMonomerSketcher(page);
      await pasteText(page, 'PEPTIDE1{A.G}$$$$V2.0');
      const before = await getExportedHelm(page);
      // Replacement intentionally clears the document before import. Select
      // Add mode here to test that a rejected molecule never enters the model.
      await clickWidget(page, 'import_btn');
      await clickPopupRow(page, 'latest', 58);
      await pasteText(page, input);
      await expectDialog(page, 'Import Error');
      expect((await widgetState(page, 'text_edit')).text).toMatch(error);
      await closeDialog(page);
      expect(await getExportedHelm(page)).toBe(before);
      await clickWidget(page, 'ala_btn');
      const atoms = await waitAtoms(page, 2);
      await page.mouse.move(atoms[0].centerX, atoms[0].centerY);
      expect(await getExportedHelm(page)).toBe(before);
    });
  }

  test('X peptides and N bases remain valid imports and survive hover and export', async ({
    page,
  }) => {
    await openMonomerSketcher(page);
    const input = 'PEPTIDE1{A.X.W}|RNA1{R(N)P}$$$$V2.0';
    await pasteText(page, input);
    const atoms = await waitAtoms(page, 6);
    await clickWidget(page, 'covalent_or_disulfide_btn');
    for (const atom of atoms.filter((item) => ['X', 'N'].includes(item.label))) {
      await page.mouse.move(atom.centerX, atom.centerY);
    }
    expect(await getExportedHelm(page)).toBe(input);
  });
});

test('SKETCH-2861: inline-SMILES HELM round trips and SMILES export preserve chemistry and chirality', async ({
  page,
  context,
}) => {
  // Jira's certification SMILES accidentally omits sulfur. Use the same
  // chemically equivalent atomistic reference as the fix's C++ regression.
  // This path never invokes monomer layout, where the corruption occurred.
  const reference = await context.newPage();
  let expectedSmiles;
  try {
    await openSketcher(reference);
    await pasteText(reference, 'O=C([C@H](CSCCN)N)O');
    expectedSmiles = await getExportedSmiles(reference);
  } finally {
    await reference.close();
  }
  const sk = await openMonomerSketcher(page);
  const input = 'PEPTIDE1{[O=C([C@H](CSCCN[H:3])N[H:1])[OH:2]]}$$$$V2.0';
  await pasteText(page, input);
  expect(await sk.copy_all_as_text('helm')).toBe(input);
  expect(await sk.copy_all_as_text('smi')).toBe(expectedSmiles);
  expect(expectedSmiles).toContain('S');
  expect(expectedSmiles).toContain('@');
  await clickWidget(page, 'cleanup_btn');
  expect(await sk.copy_all_as_text('helm')).toBe(input);
  expect(await sk.copy_all_as_text('smi')).toBe(expectedSmiles);
  await clear(page);
  await pasteText(page, input);
  expect(await getExportedHelm(page)).toBe(input);
});

test('SKETCH-2861: registered monomer symbols longer than four characters survive layout and export', async ({
  page,
}) => {
  const sk = await openMonomerSketcher(page);
  await loadDatabase(page, 'long_symbol.json');
  const input = 'PEPTIDE1{[LongAla]}$$$$V2.0';
  await pasteText(page, input);
  expect(await sk.copy_all_as_text('helm')).toBe(input);
  await clickWidget(page, 'cleanup_btn');
  expect(await sk.copy_all_as_text('helm')).toBe(input);
  expect((await scene(page)).atoms[0].label).toBe('LongAla');
  const smiles = await sk.copy_all_as_text('smi');
  expect(smiles).toContain('@');
  expect(smiles).toContain('C(=O)O');
  expect(smiles).not.toContain('*');
});
