import { expect, test } from '@playwright/test';
import {
  canvasPoint,
  clear,
  clickWidget,
  closeDialog,
  contextMenuAction,
  dialogTitle,
  expectDialog,
  getExportedHelm,
  holdPopup,
  mouseClick,
  mouseDrag,
  pasteText,
  place,
  redo,
  scene,
  screenshotView,
  selectAll,
  undo,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';
import { connect, expectAttachments, openMonomerSketcher } from './monomer_helpers.js';

test.beforeEach(async ({ page }) => {
  await openMonomerSketcher(page);
});
test.setTimeout(90_000);

test.describe('Removing bound connections — SKETCH-2865', () => {
  const DISULFIDE = 'PEPTIDE1{C}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:R3-1:R3$$$V2.0';
  const warning = (number) =>
    `${number} has been removed from this monomer but is currently bound. Continuing will remove this connection.`;

  test('editing a bound R2: Cancel returns to the edited drawing; OK removes the bond', async ({
    page,
  }) => {
    expect((await widgetState(page, 'ala_btn')).checked).toBe(true);
    await mouseDrag(page, await canvasPoint(page, -100), await canvasPoint(page, 100));
    await waitAtoms(page, 2);
    const before = await getExportedHelm(page);
    expect(before).toBe('PEPTIDE1{A.A}$$$$V2.0');
    await contextMenuAction(page, 'atom:0', 'Edit Structure...');
    await expectDialog(page, 'Edit Selected Peptide Residue');
    const r2 = (await scene(page)).atoms.find((atom) => atom.rgroup === 2);
    expect(r2).toBeTruthy();
    await clickWidget(page, 'erase_btn');
    await mouseClick(page, r2.centerX, r2.centerY);
    await expect
      .poll(async () => (await scene(page)).atoms.some((atom) => atom.rgroup === 2))
      .toBe(false);
    await closeDialog(page);
    await expectDialog(page, 'Remove Bound Connections?');
    expect((await widgetState(page, 'text_edit')).text).toBe(warning('R2'));
    expect(await getExportedHelm(page)).toBe(before);
    await closeDialog(page, 'Cancel');
    await expectDialog(page, 'Edit Selected Peptide Residue');
    expect((await scene(page)).atoms.some((atom) => atom.rgroup === 2)).toBe(false);
    await closeDialog(page);
    await expectDialog(page, 'Remove Bound Connections?');
    expect((await widgetState(page, 'text_edit')).text).toBe(warning('R2'));
    await closeDialog(page);
    expect((await waitAtoms(page, 2)).map((atom) => atom.label)).toEqual(['***', 'A']);
    const after = await getExportedHelm(page);
    expect((await scene(page)).connectors).toHaveLength(0);
    expect(after).toContain('}|PEPTIDE2{');
    await undo(page);
    expect(await getExportedHelm(page)).toBe(before);
    await redo(page);
    expect(await getExportedHelm(page)).toBe(after);
  });

  test('context mutation: Cancel retains the disulfide; OK removes it and undo restores it', async ({
    page,
  }) => {
    await clickWidget(page, 'cys_btn');
    await place(page, -100);
    await place(page, 100);
    await waitAtoms(page, 2);
    await clickWidget(page, 'covalent_or_disulfide_btn');
    await connect(page, 0, 'R3', 1, 'R3');
    expect(await getExportedHelm(page)).toBe(DISULFIDE);
    const mutate = () => contextMenuAction(page, 'atom:0', 'Mutate Residue', 'Alanine (A)', 'A');
    await mutate();
    await expectDialog(page, 'Remove Bound Connections?');
    expect((await widgetState(page, 'text_edit')).text).toBe(warning('R3'));
    await closeDialog(page, 'Cancel');
    expect(await getExportedHelm(page)).toBe(DISULFIDE);
    await mutate();
    await expectDialog(page, 'Remove Bound Connections?');
    await closeDialog(page);
    expect((await scene(page)).atoms.map((atom) => atom.label)).toEqual(['A', 'C']);
    expect((await scene(page)).connectors).toHaveLength(0);
    const after = await getExportedHelm(page);
    await undo(page);
    expect(await getExportedHelm(page)).toBe(DISULFIDE);
    await redo(page);
    expect(await getExportedHelm(page)).toBe(after);
  });

  test('compatible cysteine-to-lysine mutation preserves R3 without a warning', async ({
    page,
  }) => {
    await pasteText(page, DISULFIDE);
    await contextMenuAction(page, 'atom:0', 'Mutate Residue', 'Lysine (K)', 'K');
    await expect.poll(() => getExportedHelm(page)).toContain('PEPTIDE1{K}');
    expect(await dialogTitle(page)).toBe(null);
    expect(await getExportedHelm(page)).toContain('1:R3-1:R3');
    expect((await scene(page)).connectors).toHaveLength(1);
  });

  test('batch mutation warns once and restores both monomers and their shared link in one undo', async ({
    page,
  }) => {
    await pasteText(page, DISULFIDE);
    await selectAll(page);
    await contextMenuAction(page, 'atom:0', 'Mutate Residue', 'Alanine (A)', 'A');
    await expectDialog(page, 'Remove Bound Connections?');
    expect((await widgetState(page, 'text_edit')).text).toBe(
      'The selected replacement monomers are missing bound R3 attachment points. Continuing will remove the affected connections.',
    );
    await closeDialog(page);
    await expect.poll(() => dialogTitle(page)).toBe(null);
    expect((await scene(page)).atoms.map((atom) => atom.label)).toEqual(['A', 'A']);
    expect((await scene(page)).connectors).toHaveLength(0);
    await undo(page);
    expect(await getExportedHelm(page)).toBe(DISULFIDE);
  });
});

test.describe('Attachment sites — SKETCH-2857 / SKETCH-2834', () => {
  test('native monomers expose only their unoccupied, chemically defined sites', async ({
    page,
  }) => {
    await pasteText(page, 'PEPTIDE1{A.L.C.K.X}$$$$V2.0');
    await clickWidget(page, 'covalent_or_disulfide_btn');
    for (const [index, labels] of [
      [0, ['N']],
      [1, []],
      [2, ['S']],
      [3, ['X']],
      [4, ['C']],
    ]) {
      await expectAttachments(page, index, labels);
    }
    await screenshotView(page, 'native-attachments.png');
  });

  test('database analogs and inline SMILES derive sites and S/X labels from their structures', async ({
    page,
  }) => {
    await pasteText(
      page,
      'PEPTIDE1{[dA].[dK].[dC].[O=C([C@H](CF)N[H:1])[OH:2]].[O=C([C@H](CSCCN[H:3])N[H:1])[OH:2]]}$$$$V2.0',
    );
    await clickWidget(page, 'covalent_or_disulfide_btn');
    for (const [index, labels] of [
      [0, ['N']],
      [1, ['X']],
      [2, ['S']],
      [3, []],
      [4, ['X', 'C']],
    ]) {
      await expectAttachments(page, index, labels);
    }
  });

  test('sulfur-bearing custom monomers form a disulfide regardless of residue name', async ({
    page,
  }) => {
    await pasteText(
      page,
      'PEPTIDE1{[O=C(C(CS[H:3])N[H:1])[OH:2]]}|PEPTIDE2{[O=C(C(CS[H:3])N[H:1])[OH:2]]}$$$$V2.0',
    );
    await clickWidget(page, 'covalent_or_disulfide_btn');
    await connect(page, 0, 'R3', 1, 'R3');
    expect(await getExportedHelm(page)).toContain('1:R3-1:R3');
    await page.mouse.move(1000, 40);
    await screenshotView(page, 'custom-sulfur-disulfide.png');
  });

  test('unknown bases pair; a missing custom-nucleotide base remains usable without a crash', async ({
    page,
  }) => {
    await clickWidget(page, 'nucleic_monomer_btn');
    await clickWidget(page, 'na_custom_nt_btn');
    await place(page, -120);
    const atoms = await waitAtoms(page, 3);
    const base = atoms.find((atom) => atom.label === 'N');
    expect(base).toBeTruthy();
    await mouseClick(page, base.centerX, base.centerY);
    await waitAtoms(page, 6);
    expect(await getExportedHelm(page)).toMatch(/:pair-\d+:pair/);
    await clear(page);
    // The existing popup supports arbitrary symbols; validation belongs to a
    // future UI change. This regression guards the explicit fallback path.
    await holdPopup(page, 'na_custom_nt_btn');
    await clickWidget(page, 'base_le');
    await page.keyboard.press('ControlOrMeta+a');
    await page.keyboard.type('missing');
    await page.keyboard.press('Enter');
    await page.mouse.up();
    await place(page);
    const missing = (await waitAtoms(page, 3)).find((atom) => atom.label === 'missing');
    expect(missing).toBeTruthy();
    await clickWidget(page, 'covalent_or_disulfide_btn');
    await expectAttachments(page, missing.index, ['R2', 'H-bond']);
    expect(await getExportedHelm(page)).toContain('[missing]');
  });
});

test.describe('Phosphate backbone direction — SKETCH-2786', () => {
  for (const sugarSite of ['R1', 'R2']) {
    for (const phosphateSite of ['R1', 'R2']) {
      for (const reverse of [false, true]) {
        test(`${sugarSite} to phosphate ${phosphateSite}, ${reverse ? 'phosphate first' : 'sugar first'}`, async ({
          page,
        }) => {
          await clickWidget(page, 'nucleic_monomer_btn');
          await clickWidget(page, 'na_r_btn');
          await place(page, -100);
          await clickWidget(page, 'na_p_btn');
          await place(page, 100);
          await waitAtoms(page, 2);
          const before = await getExportedHelm(page);
          await clickWidget(page, 'covalent_or_disulfide_btn');
          await connect(
            page,
            reverse ? 1 : 0,
            reverse ? phosphateSite : sugarSite,
            reverse ? 0 : 1,
            reverse ? sugarSite : phosphateSite,
          );
          const expected = sugarSite === 'R2' ? 'RNA1{R.P}$$$$V2.0' : 'RNA1{P.R}$$$$V2.0';
          await expect.poll(() => getExportedHelm(page)).toBe(expected);
          await undo(page);
          expect(await getExportedHelm(page)).toBe(before);
          await redo(page);
          expect(await getExportedHelm(page)).toBe(expected);
        });
      }
    }
  }

  test('a partially occupied phosphate retains its existing link', async ({ page }) => {
    await pasteText(page, 'RNA1{R.P}|RNA2{R}$$$$V2.0');
    await clickWidget(page, 'covalent_or_disulfide_btn');
    await expectAttachments(page, 1, ["3'"]);
    await connect(page, 2, 'R1', 1, 'R2');
    expect(await getExportedHelm(page)).toBe('RNA1{R.P.R}$$$$V2.0');
  });
});
