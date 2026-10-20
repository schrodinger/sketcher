import { expect, test } from '@playwright/test';
import {
  canvasPoint,
  center,
  clear,
  clickPopupTool,
  clickSelector,
  clickWidget,
  closeDialog,
  contextMenuAction,
  dialogTitle,
  expectCursor,
  expectDialog,
  getExportedHelm,
  mouseClick,
  mouseDrag,
  pasteText,
  place,
  point,
  redo,
  replaceEditor,
  requireRect,
  scene,
  screenshotRegion,
  screenshotView,
  selectAll,
  undo,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';
import {
  BASE,
  BRANCHED_PEPTIDE,
  CHEM,
  PEPTIDE,
  connect,
  customMenu,
  defineCustom,
  expectAttachments,
  expectDrawingAttachments,
  openCustom,
  openMonomerSketcher,
} from './monomer_helpers.js';

test.beforeEach(async ({ page }) => {
  await openMonomerSketcher(page);
});
test.setTimeout(90_000);

test.describe('Custom menu — SKETCH-2884 / SKETCH-2886', () => {
  for (const [tab, action, title] of [
    ['amino_monomer_btn', 'Sketch Residue...', 'Define Custom Peptide Residue'],
    ['amino_monomer_btn', 'Draw Structure...', 'Define Custom Structure'],
    ['nucleic_monomer_btn', 'Sketch Base...', 'Define Custom Nucleic Acid Base'],
    ['nucleic_monomer_btn', 'Draw Structure...', 'Define Custom Structure'],
  ]) {
    test(`${tab}: ${action} opens ${title}`, async ({ page }) => {
      await clickWidget(page, tab);
      expect((await widgetState(page, 'custom_monomer_btn')).text).toBe('Custom');
      await customMenu(page);
      // Releasing the long press must leave the menu open, not invoke the
      // button's ordinary-click handler and obscure it with a dialog.
      expect(await dialogTitle(page)).toBe(null);
      expect((await scene(page)).popup.actions).toEqual([
        tab === 'amino_monomer_btn' ? 'Sketch Residue...' : 'Sketch Base...',
        '<separator>',
        'Draw Structure...',
      ]);
      await requireRect(
        page,
        `menu:${tab === 'amino_monomer_btn' ? 'Sketch Residue...' : 'Sketch Base...'}`,
      );
      await clickSelector(page, `menu:${action}`);
      await expectDialog(page, title);
      await closeDialog(page, 'Cancel');
      await expect.poll(() => dialogTitle(page)).toBe(null);
      expect((await widgetState(page, 'custom_monomer_btn')).checked).toBe(false);
      await waitAtoms(page, 0);
    });
  }

  test('ordinary click uses the current tab; button appearance follows hover and activation', async ({
    page,
  }) => {
    const rect = await requireRect(page, 'widget:custom_monomer_btn');
    await page.mouse.move(1000, 40);
    await screenshotRegion(
      page,
      'custom-palette.png',
      await requireRect(page, 'widget:monomeric_wdg'),
    );
    await screenshotRegion(page, 'custom-idle.png', rect);
    const pos = center(rect);
    await page.mouse.move(pos.x, pos.y);
    await screenshotRegion(page, 'custom-hover.png', rect);
    await clickWidget(page, 'custom_monomer_btn');
    await expectDialog(page, 'Define Custom Peptide Residue');
    await closeDialog(page, 'Cancel');
    await clickWidget(page, 'nucleic_monomer_btn');
    await clickWidget(page, 'custom_monomer_btn');
    await expectDialog(page, 'Define Custom Nucleic Acid Base');
    await replaceEditor(page, BASE);
    await closeDialog(page);
    await expect
      .poll(async () => (await widgetState(page, 'custom_monomer_btn')).checked)
      .toBe(true);
    await page.mouse.move(1000, 40);
    await screenshotRegion(
      page,
      'custom-active.png',
      await requireRect(page, 'widget:custom_monomer_btn'),
    );
    await customMenu(page);
    // This crop also records menu separators and alignment.
    await screenshotRegion(page, 'custom-base-menu.png', (await scene(page)).popup);
    await page.keyboard.press('Escape');
  });
});

test.describe('Custom structure creation', () => {
  test('SKETCH-2741: draw a custom backbone with the atomistic and R-group tools', async ({
    page,
  }) => {
    await openCustom(page);
    const left = await canvasPoint(page, -60);
    const middle = await canvasPoint(page);
    const right = await canvasPoint(page, 60);
    await clickWidget(page, 'single_bond_btn');
    await mouseDrag(page, left, middle);
    await waitAtoms(page, 2);
    await mouseDrag(page, await point(page, 1), right);
    await waitAtoms(page, 3);
    for (const [number, index] of [
      [1, 0],
      [2, 2],
    ]) {
      await clickPopupTool(page, `rgroup${number}_btn`);
      const target = await point(page, index);
      await mouseClick(page, target.x, target.y);
    }
    expect(
      (await scene(page)).atoms
        .filter((atom) => atom.rgroup)
        .map((atom) => atom.rgroup)
        .sort(),
    ).toEqual([1, 2]);
    await closeDialog(page);
    await expectCursor(page, 'custom-peptide-cursor.png');
    await place(page);
    await waitAtoms(page, 1);
    await expectDrawingAttachments(page, 0, ['N', 'C'], 'C', ['C', 'N']);
    await screenshotView(page, 'custom-peptide-attachments.png');
    const target = await point(page, 0);
    await mouseClick(page, target.x, target.y);
    await waitAtoms(page, 2);
    const helm = await getExportedHelm(page);
    expect(helm).toMatch(/^PEPTIDE1\{/);
    // CXSMILES inside a monomer can itself contain '$', so do not split HELM
    // into sections at every dollar sign.
    expect(helm).toContain('].[');
    expect(helm.endsWith('}$$$$V2.0')).toBe(true); // implicit backbone
  });

  test('SKETCH-2741: R3/R4 labels depend on the side-chain chemistry', async ({ page }) => {
    await defineCustom(page, 'peptide', '*NC(C*)(C*)C(=O)* |$_R1;;;;_R3;;_R4;;;_R2$|');
    await place(page);
    await waitAtoms(page, 1);
    await expectDrawingAttachments(page, 0, ['N', 'C', 'X', 'R4'], 'C', ['C', 'N']);
    await clear(page);
    await defineCustom(page, 'peptide', '*NC(CS*)C(=O)* |$_R1;;;;;_R3;;;_R2$|');
    await place(page);
    await waitAtoms(page, 1);
    await expectDrawingAttachments(page, 0, ['N', 'C', 'S'], 'C', ['C', 'N']);
  });

  test('SKETCH-2825: custom bases expose a glycosidic site and create a hydrogen pair', async ({
    page,
  }) => {
    await defineCustom(page, 'base', BASE);
    await expectCursor(page, 'custom-base-cursor.png');
    await place(page);
    await waitAtoms(page, 1);
    await expectDrawingAttachments(page, 0, ['N1/9', 'H-bond'], 'H-bond', ['H-bond']);
    await screenshotView(page, 'custom-base-attachments.png');
    const target = await point(page, 0);
    await mouseClick(page, target.x, target.y);
    await waitAtoms(page, 2);
    expect(await getExportedHelm(page)).toMatch(/1:pair-1:pair/);
  });

  for (const tab of ['amino_monomer_btn', 'nucleic_monomer_btn']) {
    test(`SKETCH-2771: define and connect CHEM from ${tab}`, async ({ page }) => {
      await clickWidget(page, tab);
      await defineCustom(page, 'chem', CHEM);
      await expectCursor(page, `custom-chem-cursor-${tab}.png`);
      await place(page, -100);
      await place(page, 100);
      await waitAtoms(page, 2);
      await expectDrawingAttachments(page, 0, ['R1', 'R2', 'R3'], 'R1', ['R1', 'R1']);
      await screenshotView(page, `custom-chem-${tab}.png`);
      await clickWidget(page, 'covalent_or_disulfide_btn');
      await connect(page, 0, 'R2', 1, 'R1');
      const helm = await getExportedHelm(page);
      expect(helm).toMatch(/^CHEM1\{\[.*\]\}\|CHEM2\{\[.*\]\}\$/);
      expect(helm).toMatch(/CHEM1,CHEM2,1:R2-1:R1/);
    });
  }

  test('SKETCH-2741: accepting a custom residue mutates the current selection', async ({
    page,
  }) => {
    await pasteText(page, 'PEPTIDE1{A.G}$$$$V2.0');
    await selectAll(page);
    await defineCustom(page, 'peptide', PEPTIDE);
    const atoms = await waitAtoms(page, 2);
    expect(atoms.map((atom) => atom.label)).toEqual(['***', '***']);
    expect(await getExportedHelm(page)).toContain('].[');
    await undo(page);
    expect(await getExportedHelm(page)).toBe('PEPTIDE1{A.G}$$$$V2.0');
  });

  test('SKETCH-2826: empty/disconnected drawings cannot be accepted; duplicate R groups explain the error', async ({
    page,
  }) => {
    await openCustom(page);
    expect((await requireRect(page, 'dialog:OK')).enabled).toBe(false);
    await replaceEditor(page, 'CC.C');
    expect((await requireRect(page, 'dialog:OK')).enabled).toBe(false);
    await replaceEditor(page, '*C* |$_R1;;_R1$|');
    await closeDialog(page);
    await expectDialog(page, 'Invalid Attachment Points');
    expect((await widgetState(page, 'text_edit')).text).toBe(
      'Multiple R1 attachment points found. All attachment points must be unique.',
    );
    await closeDialog(page);
    await expectDialog(page, 'Define Custom Peptide Residue');
    await replaceEditor(page, PEPTIDE);
    await closeDialog(page);
    await place(page);
    await waitAtoms(page, 1);
  });
});

test.describe('Existing structure editing — SKETCH-2826 / SKETCH-2886', () => {
  for (const [type, smiles, title] of [
    ['peptide', PEPTIDE, 'Edit Selected Peptide Residue'],
    ['base', BASE, 'Edit Selected Nucleic Acid Base'],
    ['chem', CHEM, 'Edit Selected Structure'],
  ]) {
    test(`${type}: editor reopens the saved structure and Cancel discards atom edits`, async ({
      page,
    }) => {
      await defineCustom(page, type, smiles);
      await place(page);
      await waitAtoms(page, 1);
      const before = await getExportedHelm(page);
      await contextMenuAction(page, 'atom:0', 'Edit Structure...');
      await expectDialog(page, title);
      const original = (await scene(page)).atoms.map(({ element, rgroup }) => ({
        element,
        rgroup,
      }));
      expect(original.length).toBeGreaterThan(2);
      await replaceEditor(page, type === 'peptide' ? BRANCHED_PEPTIDE : '*CCC* |$_R1;;;;_R2$|');
      await closeDialog(page, 'Cancel');
      expect(await getExportedHelm(page)).toBe(before);
      await contextMenuAction(page, 'atom:0', 'Edit Structure...');
      await expectDialog(page, title);
      expect((await scene(page)).atoms.map(({ element, rgroup }) => ({ element, rgroup }))).toEqual(
        original,
      );
      await replaceEditor(page, type === 'peptide' ? BRANCHED_PEPTIDE : '*CCC* |$_R1;;;;_R2$|');
      await closeDialog(page);
      const after = await getExportedHelm(page);
      expect(after).not.toBe(before);
      if (type === 'peptide') {
        await clickWidget(page, 'covalent_or_disulfide_btn');
        await expectAttachments(page, 0, ['N', 'C', 'X']);
      }
      await undo(page);
      expect(await getExportedHelm(page)).toBe(before);
      await redo(page);
      expect(await getExportedHelm(page)).toBe(after);
    });
  }

  test('native alanine editing preserves the backbone and becomes a custom bead', async ({
    page,
  }) => {
    await pasteText(page, 'PEPTIDE1{A.G}$$$$V2.0');
    await contextMenuAction(page, 'atom:0', 'Edit Structure...');
    await expectDialog(page, 'Edit Selected Peptide Residue');
    const atoms = (await scene(page)).atoms;
    expect(
      atoms
        .filter((atom) => atom.rgroup)
        .map((atom) => atom.rgroup)
        .sort(),
    ).toEqual([1, 2]);
    expect(atoms.some((atom) => atom.element === 'N')).toBe(true);
    await replaceEditor(page, BRANCHED_PEPTIDE);
    await closeDialog(page);
    expect((await waitAtoms(page, 2)).map((atom) => atom.label)).toEqual(['***', 'G']);
    expect(await getExportedHelm(page)).toMatch(/\]\.G\}/);
  });
});
