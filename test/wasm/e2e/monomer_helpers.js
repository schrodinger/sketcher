import { expect } from '@playwright/test';
import path from 'node:path';
import {
  FIXTURES,
  clickSelector,
  clickWidget,
  closeDialog,
  dialogTitle,
  expectDialog,
  getRect,
  holdPopup,
  hoverAtom,
  mouseDrag,
  openSketcher,
  replaceEditor,
  scene,
  widgetState,
} from './e2e_helpers.js';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { setFilePickerResult } from './sketcher/wrappers/sketcher_wasm.js';

export const PEPTIDE = '*NCC(=O)* |$_R1;;;;;_R2$|';
export const BRANCHED_PEPTIDE = '*NC(C*)C(=O)* |$_R1;;;;_R3;;;_R2$|';
export const BASE = '*n1cnc2c(N)ncnc21 |$_R1;;;;;;;;;;$|';
export const CHEM = '*CC(*)C* |$_R1;;;_R2;;_R3$|';

export async function openMonomerSketcher(page) {
  const sk = await openSketcher(page);
  await clickWidget(page, 'monomeric_btn');
  return sk;
}

export async function customMenu(page) {
  await holdPopup(page, 'custom_monomer_btn');
  await page.mouse.up();
  await expect.poll(async () => Boolean(await getRect(page, 'menu:Draw Structure...'))).toBe(true);
}

export async function openCustom(page, type = 'peptide') {
  if (type !== 'chem') {
    await clickWidget(page, type === 'peptide' ? 'amino_monomer_btn' : 'nucleic_monomer_btn');
  }
  await customMenu(page);
  await clickSelector(
    page,
    `menu:${type === 'chem' ? 'Draw Structure...' : type === 'peptide' ? 'Sketch Residue...' : 'Sketch Base...'}`,
  );
  await expectDialog(
    page,
    `Define Custom ${type === 'peptide' ? 'Peptide Residue' : type === 'base' ? 'Nucleic Acid Base' : 'Structure'}`,
  );
}

export async function defineCustom(page, type, smiles) {
  await openCustom(page, type);
  await replaceEditor(page, smiles);
  await closeDialog(page);
  await expect.poll(() => dialogTitle(page)).toBe(null);
}

export async function attachments(page, index) {
  await hoverAtom(page, index);
  // Hover items are created on a Qt event turn after mouse movement.
  await page.waitForTimeout(100);
  return (await scene(page)).attachments.filter((item) => item.owner === index);
}

export async function expectAttachments(page, index, labels) {
  await expect
    .poll(async () => (await attachments(page, index)).map((item) => item.text).sort())
    .toEqual([...labels].sort());
}

export async function expectDrawingAttachments(page, index, labels, previewLabel, hintLabels) {
  // Certification counts the default site's label on the hint connector.
  // Check the actual painted hint text as well as every remaining nubbin;
  // hydrogen bonds use one shared label instead of two endpoint labels.
  expect(labels).toContain(previewLabel);
  expect(hintLabels).toContain(previewLabel);
  await hoverAtom(page, index);
  await expect
    .poll(async () => {
      const state = await scene(page);
      expect(
        state,
        'WASM artifact must include hintConnections in the scene bridge',
      ).toHaveProperty('hintConnections');
      return {
        unbound: state.attachments
          .filter((item) => item.owner === index)
          .map((item) => item.text)
          .sort(),
        hints: state.hintConnections
          .filter((hint) => hint.owners.includes(index))
          .map((hint) => ({
            owners: [...hint.owners].sort((a, b) => a - b),
            labels: hint.labels.map((item) => item.text).sort(),
          })),
      };
    })
    .toEqual({
      unbound: labels.filter((label) => label !== previewLabel).sort(),
      hints: [{ owners: [index], labels: [...hintLabels].sort() }],
    });
}

export async function attachment(page, index, name) {
  let result;
  await expect
    .poll(
      async () => {
        result = (await attachments(page, index)).find((item) => item.name === name);
        return Boolean(result);
      },
      { message: `atom ${index} exposes ${name}` },
    )
    .toBe(true);
  return { x: result.hitX, y: result.hitY };
}

export async function connect(page, first, firstName, second, secondName) {
  const start = await attachment(page, first, firstName);
  const end = await attachment(page, second, secondName);
  await mouseDrag(page, start, end);
}

export async function loadDatabase(page, filename) {
  await setFilePickerResult(page, path.join(FIXTURES, filename));
  await new Sketcher(page).import_menu('load_monomer_database');
  await expectDialog(page, 'Monomer Database', 30000);
  await expect
    .poll(async () => (await widgetState(page, 'text_edit')).text)
    .toBe(`${filename} read successfully.`);
  await closeDialog(page);
  await expect.poll(() => dialogTitle(page)).toBe(null);
}
