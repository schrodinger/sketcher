import { expect, test } from '@playwright/test';
import path from 'node:path';
import { setFilePickerResult } from './sketcher/wrappers/sketcher_wasm.js';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import {
  FIXTURES,
  choosePopup,
  clear,
  clickSelector,
  clickWidget,
  closeDialog,
  dismissPopup,
  expectCursor,
  expectDialog,
  getExportedHelm,
  getRect,
  isWidgetVisible,
  pasteText,
  place,
  popupButtons,
  popupInfo,
  requireRect,
  screenshotView,
  selectAll,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';
import {
  MONOMER_DB_TABLE,
  chooseFromMonomerDatabase,
  closeMonomerDatabase,
  expectAttachments,
  loadDatabase,
  findMonomerDatabaseRow,
  monomerDatabaseRow,
  monomerDatabaseSymbols,
  monomerDatabaseTable,
  openMonomerDatabase,
  openMonomerSketcher,
  searchMonomerDatabase,
} from './monomer_helpers.js';

test.beforeEach(async ({ page }) => {
  await openMonomerSketcher(page);
});
test.setTimeout(120_000);

async function expectPalette(page, button, symbols) {
  expect((await popupButtons(page, button)).map((item) => item.text)).toEqual(symbols);
  await dismissPopup(page);
}

/**
 * Check the symbols in the Monomer Database popup's table on the current page,
 * with no search text
 */
async function expectTable(page, symbols) {
  await openMonomerDatabase(page);
  await searchMonomerDatabase(page, '');
  expect(await monomerDatabaseSymbols(page)).toEqual(symbols);
  await closeMonomerDatabase(page);
}

test('SKETCH-2870: defaults and successive database replacements refresh the Monomer Database table', async ({
  page,
}) => {
  // with no custom database, the table holds the core monomers with no
  // natural analog
  await expectTable(page, ['ac', 'am']);
  await chooseFromMonomerDatabase(page, 'ac');
  await expectCursor(page, 'monomer-db-ac-cursor.png');
  await place(page);
  expect((await waitAtoms(page, 1))[0].label).toBe('ac');
  await clickWidget(page, 'ala_btn');
  await expectAttachments(page, 0, ['C']);
  await screenshotView(page, 'monomer-db-ac.png');
  await clear(page);

  // the core database has no nucleic acids or structures without an analog
  await clickWidget(page, 'nucleic_monomer_btn');
  await expectTable(page, []);
  await loadDatabase(page, 'chembl_monomers_no_id.json');
  await clickWidget(page, 'amino_monomer_btn');
  await openMonomerDatabase(page);
  await searchMonomerDatabase(page, '');
  expect((await monomerDatabaseTable(page)).rowCount).toBeGreaterThan(100);
  // the core monomers with no natural analog are still in the table
  await searchMonomerDatabase(page, 'ac');
  expect(await findMonomerDatabaseRow(page, 'ac')).toBeGreaterThanOrEqual(0);
  await closeMonomerDatabase(page);

  await loadDatabase(page, 'chem_monomers.json');
  await expectTable(page, ['ac', 'am', 'CHM1', 'CHM2', 'CHM3']);
  await chooseFromMonomerDatabase(page, 'CHM2');
  await expectCursor(page, 'monomer-db-chem-cursor.png');
  await place(page, -100);
  expect((await waitAtoms(page, 1))[0].label).toBe('CHM2');
  expect(await getExportedHelm(page)).toContain('CHEM1{[CHM2]}');
  await screenshotView(page, 'monomer-db-chem.png');
  await clear(page);
  await clickWidget(page, 'nucleic_monomer_btn');
  await expectTable(page, ['CHM1', 'CHM2', 'CHM3']);
  await chooseFromMonomerDatabase(page, 'CHM1');
  await expectCursor(page, 'nucleic-monomer-db-chem-cursor.png');
  await place(page);
  expect(await getExportedHelm(page)).toContain('CHEM1{[CHM1]}');
  await clear(page);

  await loadDatabase(page, 'na_monomers.json');
  await expectTable(page, ['NA1', 'NA2']);
  await chooseFromMonomerDatabase(page, 'NA1');
  await expectCursor(page, 'monomer-db-base-cursor.png');
  await place(page);
  expect((await waitAtoms(page, 1))[0].label).toBe('NA1');
  await screenshotView(page, 'monomer-db-base.png');
  await clear(page);
  await loadDatabase(page, 'empty_monomers.json');
  await expectTable(page, []);
  await clickWidget(page, 'amino_monomer_btn');
  await expectTable(page, ['ac', 'am']);
});

test('SKETCH-2870: scroll the Monomer Database table to a real off-screen monomer and draw it', async ({
  page,
}) => {
  await loadDatabase(page, 'chembl_monomers_no_id.json');
  await openMonomerDatabase(page);
  const table = await monomerDatabaseTable(page);
  const lastRow = table.rowCount - 1;
  expect(lastRow).toBeGreaterThan(100);
  // the bridge scrolls the row into view before reporting where it is
  const last = await monomerDatabaseRow(page, lastRow);
  expect(last.y).toBeGreaterThanOrEqual(table.y);
  expect(last.y + last.height).toBeLessThanOrEqual(table.y + table.height);
  const symbol = last.cells[0];
  await clickSelector(page, `table_row:${MONOMER_DB_TABLE}:${lastRow}`);
  await expect
    .poll(async () => (await widgetState(page, 'monomer_db_last_used_btn')).text)
    .toBe(symbol);
  await closeMonomerDatabase(page);
  await place(page);
  expect((await waitAtoms(page, 1))[0].label).toBeTruthy();
  expect(await getExportedHelm(page)).toContain(symbol);
});

test('SKETCH-2870: each tab remembers its own database monomer; tab changes preserve selected content', async ({
  page,
}) => {
  await loadDatabase(page, 'chem_monomers.json');
  await chooseFromMonomerDatabase(page, 'CHM2');
  await clickWidget(page, 'nucleic_monomer_btn');
  await chooseFromMonomerDatabase(page, 'CHM1');
  await clickWidget(page, 'amino_monomer_btn');
  await expect
    .poll(async () => (await widgetState(page, 'monomer_db_last_used_btn')).text)
    .toBe('CHM2');
  await place(page, -120);
  expect((await waitAtoms(page, 1))[0].label).toBe('CHM2');
  await clickWidget(page, 'nucleic_monomer_btn');
  await expect
    .poll(async () => (await widgetState(page, 'monomer_db_last_used_btn')).text)
    .toBe('CHM1');
  await place(page, 120);
  expect((await waitAtoms(page, 2)).map((atom) => atom.label)).toEqual(['CHM2', 'CHM1']);
  const before = await getExportedHelm(page);
  await selectAll(page);
  await clickWidget(page, 'amino_monomer_btn');
  await clickWidget(page, 'nucleic_monomer_btn');
  expect(await getExportedHelm(page)).toBe(before);
});

test('SKETCH-2870: choosing a database peptide mutates the selection rather than placing a bead', async ({
  page,
}) => {
  await pasteText(page, 'PEPTIDE1{A}$$$$V2.0');
  await selectAll(page);
  await chooseFromMonomerDatabase(page, 'ac');
  expect((await waitAtoms(page, 1))[0].label).toBe('ac');
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{[ac]}$$$$V2.0');
});

test('SKETCH-2788: the upload menu reports success and replaces large analog palettes with native definitions', async ({
  page,
}) => {
  await loadDatabase(page, 'chembl_monomers_no_id.json');
  const large = await popupButtons(page, 'gly_btn');
  expect(large.length).toBeGreaterThan(3);
  expect(large.some((button) => button.text === 'Boc_G')).toBe(true);
  await dismissPopup(page);
  await choosePopup(page, 'gly_btn', 'Boc_G', { scroll: true });
  await place(page);
  expect(await getExportedHelm(page)).toContain('{[Boc_G]}');
  await clear(page);
  await loadDatabase(page, 'native_monomer_definitions.json');
  expect((await widgetState(page, 'gly_btn')).text).toBe('G');
  await expectPalette(page, 'gly_btn', ['G', 'meG']);
  await choosePopup(page, 'gly_btn', 'meG');
  await place(page);
  expect(await getExportedHelm(page)).toContain('{[meG]}');
});

test('SKETCH-2788: malformed JSON explains the failure; picker cancellation leaves content and palettes intact', async ({
  page,
}) => {
  await pasteText(page, 'PEPTIDE1{A}$$$$V2.0');
  const before = await getExportedHelm(page);
  await setFilePickerResult(page, path.join(FIXTURES, 'invalid_monomers.json'));
  await new Sketcher(page).import_menu('load_monomer_database');
  await expectDialog(page, 'Monomer Database Error');
  expect((await widgetState(page, 'text_edit')).text).toMatch(/parse|JSON|invalid/i);
  await closeDialog(page);
  expect(await getExportedHelm(page)).toBe(before);
  await expectTable(page, ['ac', 'am']);
  await page.evaluate(() => {
    window.showOpenFilePicker = async () => {
      throw new DOMException('Picker cancelled', 'AbortError');
    };
  });
  await new Sketcher(page).import_menu('load_monomer_database');
  // Reopening the popup proves the asynchronous cancellation returned control
  // to the UI. No fixed delay or assertion of silence alone is used.
  await expectTable(page, ['ac', 'am']);
  expect(await getExportedHelm(page)).toBe(before);
});

test('SKETCH-2827: a mutation from the Monomer Database popup is remembered as the DB monomer of its page without changing the page tool', async ({
  page,
}) => {
  await pasteText(page, 'PEPTIDE1{A}$$$$V2.0');
  await selectAll(page);
  await chooseFromMonomerDatabase(page, 'ac');
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{[ac]}$$$$V2.0');
  // Switching pages and back with the selection still active changes nothing
  await clickWidget(page, 'nucleic_monomer_btn');
  await clickWidget(page, 'amino_monomer_btn');
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{[ac]}$$$$V2.0');
  await clickWidget(page, 'clear_selection_btn');

  // The mutation didn't change the page's last monomeric tool. Clearing the
  // selection leaves the select tool active, and a page switch keeps a shared
  // tool like that, so choose a monomer on the NA page before switching back.
  await clickWidget(page, 'nucleic_monomer_btn');
  await clickWidget(page, 'na_a_btn');
  await clickWidget(page, 'amino_monomer_btn');
  await expect.poll(async () => (await widgetState(page, 'ala_btn')).checked).toBe(true);

  // The last-used button reactivates the remembered DB monomer directly
  expect((await widgetState(page, 'monomer_db_last_used_btn')).text).toBe('ac');
  await clickWidget(page, 'monomer_db_last_used_btn');
  await expect
    .poll(async () => (await widgetState(page, 'monomer_db_last_used_btn')).checked)
    .toBe(true);
  expect(await isWidgetVisible(page, 'monomer_db_popup')).toBe(false);
  await place(page, 0, 120);
  expect((await waitAtoms(page, 2)).map((atom) => atom.label)).toEqual(['ac', 'ac']);
});

test('SKETCH-2827: loading a database without the active DB monomer falls back to the standard tool of the page', async ({
  page,
}) => {
  await loadDatabase(page, 'chem_monomers.json');
  await chooseFromMonomerDatabase(page, 'CHM2');
  expect((await widgetState(page, 'monomer_db_last_used_btn')).checked).toBe(true);
  // na_monomers.json has no CHEM monomers
  await loadDatabase(page, 'na_monomers.json');
  await expect.poll(async () => (await widgetState(page, 'ala_btn')).checked).toBe(true);
  expect((await widgetState(page, 'monomer_db_last_used_btn')).checked).toBe(false);
  expect((await widgetState(page, 'monomer_db_last_used_btn')).text).toBe('');
  await place(page);
  expect((await waitAtoms(page, 1))[0].label).toBe('A');
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{A}$$$$V2.0');
  await clear(page);

  await clickWidget(page, 'nucleic_monomer_btn');
  await chooseFromMonomerDatabase(page, 'NA1');
  expect((await widgetState(page, 'monomer_db_last_used_btn')).checked).toBe(true);
  // chem_monomers.json has no RNA monomers
  await loadDatabase(page, 'chem_monomers.json');
  await expect.poll(async () => (await widgetState(page, 'na_rna_btn')).checked).toBe(true);
  expect((await widgetState(page, 'monomer_db_last_used_btn')).checked).toBe(false);
  await place(page);
  await waitAtoms(page, 3);
  expect(await getExportedHelm(page)).toBe('RNA1{R(A)P}$$$$V2.0');
});

test('SKETCH-2827: a long-press popup with no room beside its button in a small canvas is clamped inside it', async ({
  page,
}) => {
  // gly_btn is in the side bar's right-hand column, and its ChEMBL palette is
  // one of the largest, so a small canvas can still hold the palette while
  // leaving it no room to the right of the button's left edge or on either
  // side of the button vertically
  await loadDatabase(page, 'chembl_monomers_no_id.json');
  const full = await popupInfo(page, 'gly_btn');
  await dismissPopup(page);
  const width = full.width + 26;
  const height = full.height + 76;
  const viewBefore = await requireRect(page, 'widget:view');
  await page.setViewportSize({ width, height });
  // Qt only sees the new canvas size once the browser's ResizeObserver fires
  await expect
    .poll(async () => {
      const view = await getRect(page, 'widget:view');
      return (
        view !== null && (view.width !== viewBefore.width || view.height !== viewBefore.height)
      );
    })
    .toBe(true);

  // Check that the shrunken canvas forces both kinds of clamping, using the
  // button's position after the resize
  const button = await requireRect(page, 'widget:gly_btn');
  const roomAbove = button.y;
  const roomBelow = height - (button.y + button.height);
  expect(button.x + full.width).toBeGreaterThan(width);
  expect(roomAbove).toBeLessThan(full.height);
  expect(roomBelow).toBeLessThan(full.height);

  const popup = await popupInfo(page, 'gly_btn');
  // A long-press popup keeps its size, opens on the side of the button with
  // more room, and is then clamped into the canvas
  expect(popup.width).toBe(full.width);
  expect(popup.height).toBe(full.height);
  expect(popup.x).toBe(width - full.width);
  expect(popup.y).toBe(roomBelow >= roomAbove ? height - full.height : 0);
  await dismissPopup(page);
});
