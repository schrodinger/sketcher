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
  hoverAtom,
  isWidgetVisible,
  mouseDrag,
  openSketcher,
  replaceEditor,
  requireRect,
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
  await clickWidget(page, 'custom_monomer_btn');
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

/** The objectName of the Monomer Database popup's table */
export const MONOMER_DB_TABLE = 'monomer_db_table';

export async function isMonomerDatabaseOpen(page) {
  return isWidgetVisible(page, 'monomer_db_popup');
}

/**
 * Open the Monomer Database popup from its side bar button, unless it's
 * already open
 */
export async function openMonomerDatabase(page) {
  if (await isMonomerDatabaseOpen(page)) return;
  await clickWidget(page, 'monomer_db_btn');
  await expect.poll(() => isMonomerDatabaseOpen(page)).toBe(true);
}

/**
 * Close the Monomer Database popup with its Close link, unless it's already
 * closed. A click outside the popup only dismisses it, so a test that needs to
 * click the canvas next should close the popup first.
 */
export async function closeMonomerDatabase(page) {
  if (!(await isMonomerDatabaseOpen(page))) return;
  await clickWidget(page, 'monomer_db_close_link');
  await expect.poll(() => isMonomerDatabaseOpen(page)).toBe(false);
}

/** Replace the open popup's search text by typing */
export async function searchMonomerDatabase(page, text) {
  await clickWidget(page, 'monomer_db_search_le');
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.press('Backspace');
  if (text) await page.keyboard.type(text, { delay: 10 });
  await expect.poll(async () => (await widgetState(page, 'monomer_db_search_le')).text).toBe(text);
}

/**
 * Return the open popup's table state: "rowCount", "headers", "selectedRow" and
 * "currentRow"
 */
export async function monomerDatabaseTable(page) {
  return requireRect(page, `widget:${MONOMER_DB_TABLE}`);
}

/**
 * Return a row of the open popup's table, after scrolling it into view, with
 * its "cells" (the visible columns' text), "selected" and "current" state
 */
export async function monomerDatabaseRow(page, row) {
  return requireRect(page, `table_row:${MONOMER_DB_TABLE}:${row}`);
}

/** Return the symbols in the open popup's table, in display order */
export async function monomerDatabaseSymbols(page) {
  const { rowCount } = await monomerDatabaseTable(page);
  const symbols = [];
  for (let row = 0; row < rowCount; row += 1) {
    symbols.push((await monomerDatabaseRow(page, row)).cells[0]);
  }
  return symbols;
}

/**
 * Return the row of the open popup's table that shows the given symbol, and
 * the given Type text ("Peptide", "Nucleic acid" or "Structure") if any, or -1
 */
export async function findMonomerDatabaseRow(page, symbol, type) {
  const { rowCount } = await monomerDatabaseTable(page);
  for (let row = 0; row < rowCount; row += 1) {
    const { cells } = await monomerDatabaseRow(page, row);
    if (cells[0] === symbol && (type === undefined || cells.includes(type))) {
      return row;
    }
  }
  return -1;
}

/**
 * Choose a monomer from the Monomer Database popup by searching for it and
 * clicking its row. Without a canvas selection this activates the monomer, and
 * with one it mutates the selection.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} symbol - the monomer's symbol
 * @param {object} [options]
 * @param {string} [options.search] - the search text, the symbol by default
 * @param {string} [options.type] - the Type column text, to tell apart monomers
 *   that share a symbol
 * @param {boolean} [options.close] - whether to close the popup afterwards
 */
export async function chooseFromMonomerDatabase(
  page,
  symbol,
  { search = symbol, type, close = true } = {},
) {
  await openMonomerDatabase(page);
  await searchMonomerDatabase(page, search);
  let row = -1;
  await expect
    .poll(
      async () => {
        row = await findMonomerDatabaseRow(page, symbol, type);
        return row;
      },
      { message: `the Monomer Database lists ${symbol}` },
    )
    .toBeGreaterThanOrEqual(0);
  await clickSelector(page, `table_row:${MONOMER_DB_TABLE}:${row}`);
  if (close) await closeMonomerDatabase(page);
}
