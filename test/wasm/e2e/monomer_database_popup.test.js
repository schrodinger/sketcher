import { expect, test } from '@playwright/test';
import path from 'node:path';
import { setFilePickerResult } from './sketcher/wrappers/sketcher_wasm.js';
import {
  FIXTURES,
  clear,
  clickSelector,
  clickWidget,
  closeDialog,
  dialogTitle,
  expectDialog,
  getExportedHelm,
  getRect,
  mouseClick,
  pasteText,
  place,
  requireRect,
  screenshotRegion,
  selectAll,
  undo,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';
import {
  MONOMER_DB_TABLE,
  chooseFromMonomerDatabase,
  closeMonomerDatabase,
  findMonomerDatabaseRow,
  isMonomerDatabaseOpen,
  loadDatabase,
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

async function matchCount(page) {
  return (await widgetState(page, 'monomer_db_match_count_lbl')).text;
}

async function lastUsed(page) {
  return widgetState(page, 'monomer_db_last_used_btn');
}

async function searchText(page) {
  return (await widgetState(page, 'monomer_db_search_le')).text;
}

/** Click the link that loads or updates the monomer database */
async function clickLoadLink(page, expectedText) {
  expect((await widgetState(page, 'monomer_db_load_link')).text).toBe(expectedText);
  await clickWidget(page, 'monomer_db_load_link');
  // The popup closes before the file chooser opens
  await expect.poll(() => isMonomerDatabaseOpen(page)).toBe(false);
}

test('SKETCH-2827: the popup opens from its button and closes with the x, Close, its button and an outside click', async ({
  page,
}) => {
  const closers = {
    x: () => clickWidget(page, 'monomer_db_close_btn'),
    close: () => clickWidget(page, 'monomer_db_close_link'),
    button: () => clickWidget(page, 'monomer_db_btn'),
    outside: async () => {
      const view = await requireRect(page, 'widget:view');
      await mouseClick(page, view.x + view.width / 2, view.y - 10);
    },
  };
  for (const [name, close] of Object.entries(closers)) {
    await openMonomerDatabase(page);
    await close();
    await expect
      .poll(() => isMonomerDatabaseOpen(page), { message: `closed with ${name}` })
      .toBe(false);
  }
  // Opening the popup never changes the tool
  expect((await widgetState(page, 'ala_btn')).checked).toBe(true);
  expect((await lastUsed(page)).checked).toBe(false);
});

test('SKETCH-2827 / SKETCH-2929: a canvas click outside the popup dismisses it without placing a monomer', async ({
  page,
}) => {
  await chooseFromMonomerDatabase(page, 'ac', { close: false });
  expect((await lastUsed(page)).checked).toBe(true);
  // A click that dismisses a side bar popup only dismisses it, on every
  // platform
  await place(page);
  await expect.poll(() => isMonomerDatabaseOpen(page)).toBe(false);
  await waitAtoms(page, 0);
  // Click somewhere else, so that this isn't treated as a double-click
  await place(page, 80);
  expect((await waitAtoms(page, 1))[0].label).toBe('ac');
});

test('SKETCH-2827: with no database the table shows ac and am, and the Load Database link loads one', async ({
  page,
}) => {
  await openMonomerDatabase(page);
  expect(await monomerDatabaseSymbols(page)).toEqual(['ac', 'am']);
  expect(await matchCount(page)).toBe('2 matches');
  expect((await monomerDatabaseTable(page)).headers).toEqual(['Symbol', 'Name', 'Type', 'Analog']);

  await setFilePickerResult(page, path.join(FIXTURES, 'mixed_monomers.json'));
  await clickLoadLink(page, 'Load Database...');
  await expectDialog(page, 'Monomer Database', 30000);
  await expect
    .poll(async () => (await widgetState(page, 'text_edit')).text)
    .toBe('mixed_monomers.json read successfully.');
  await closeDialog(page);
  await expect.poll(() => dialogTitle(page)).toBe(null);
  expect(await isMonomerDatabaseOpen(page)).toBe(false);

  await openMonomerDatabase(page);
  expect((await widgetState(page, 'monomer_db_load_link')).text).toBe('Update Database...');
  // the AA page's default filter shows peptides and structures
  expect(await monomerDatabaseSymbols(page)).toEqual([
    'ac',
    'Aib',
    'am',
    'Bz',
    'CHM1',
    'CHM2',
    'dAla',
    'Hyp',
    'LongPeptideSymbol',
  ]);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: loading or updating closes the popup, which stays closed and keeps its search and sort', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabase(page);
  await searchMonomerDatabase(page, 'a');
  // sort by name, descending
  const table = await monomerDatabaseTable(page);
  const nameHeader = {
    x: table.x + table.width / 2,
    y: table.y + 10,
  };
  await mouseClick(page, nameHeader.x, nameHeader.y);
  await mouseClick(page, nameHeader.x, nameHeader.y);
  const sorted = await monomerDatabaseSymbols(page);

  // cancelling the file chooser
  await page.evaluate(() => {
    window.showOpenFilePicker = async () => {
      throw new DOMException('Picker cancelled', 'AbortError');
    };
  });
  await clickLoadLink(page, 'Update Database...');
  // Reopening proves that the cancellation returned control to the UI
  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('a');
  expect(await monomerDatabaseSymbols(page)).toEqual(sorted);

  // a file that can't be read
  await setFilePickerResult(page, path.join(FIXTURES, 'invalid_monomers.json'));
  await clickLoadLink(page, 'Update Database...');
  await expectDialog(page, 'Monomer Database Error');
  await closeDialog(page);
  expect(await isMonomerDatabaseOpen(page)).toBe(false);

  // a successful load shows the new rows, with the same search and sort
  await setFilePickerResult(page, path.join(FIXTURES, 'chem_monomers.json'));
  await openMonomerDatabase(page);
  await clickLoadLink(page, 'Update Database...');
  await expectDialog(page, 'Monomer Database', 30000);
  await closeDialog(page);
  expect(await isMonomerDatabaseOpen(page)).toBe(false);
  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('a');
  // only the core monomers have an "a" in their symbol or name
  expect((await monomerDatabaseSymbols(page)).sort()).toEqual(['ac', 'am']);
  expect((await monomerDatabaseTable(page)).headers).toEqual(['Symbol', 'Name', 'Type', 'Analog']);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: search narrows the table and updates the match count', async ({ page }) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabase(page);
  expect(await matchCount(page)).toBe('9 matches');
  // a substring of a name, case-insensitively
  await searchMonomerDatabase(page, 'HYDROXY');
  expect(await monomerDatabaseSymbols(page)).toEqual(['Hyp']);
  expect(await matchCount(page)).toBe('1 match');
  // a substring of a symbol
  await searchMonomerDatabase(page, 'chm');
  expect(await monomerDatabaseSymbols(page)).toEqual(['CHM1', 'CHM2']);
  expect(await matchCount(page)).toBe('2 matches');
  await searchMonomerDatabase(page, 'no such monomer');
  expect((await monomerDatabaseTable(page)).rowCount).toBe(0);
  expect(await matchCount(page)).toBe('0 matches');
  await screenshotRegion(page, 'monomer-db-no-matches.png', await monomerDatabaseTable(page));
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: a chosen monomer is drawn, labels the last-used button and is reactivated by it', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await chooseFromMonomerDatabase(page, 'dAla', { close: false });
  // choosing a row keeps the popup open and selects the row
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  const table = await monomerDatabaseTable(page);
  expect(table.selectedRow).toBe(await findMonomerDatabaseRow(page, 'dAla'));
  await closeMonomerDatabase(page);
  const button = await lastUsed(page);
  expect(button.text).toBe('dAla');
  expect(button.toolTip).toBe('D-Alanine');
  expect(button.checked).toBe(true);
  expect((await widgetState(page, 'ala_btn')).checked).toBe(false);
  await place(page);
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{[dAla]}$$$$V2.0');
  await clear(page);

  // another tool, then the last-used button
  await clickWidget(page, 'gly_btn');
  expect((await lastUsed(page)).checked).toBe(false);
  await clickWidget(page, 'monomer_db_last_used_btn');
  await expect.poll(async () => (await lastUsed(page)).checked).toBe(true);
  expect(await isMonomerDatabaseOpen(page)).toBe(false);
  await place(page);
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{[dAla]}$$$$V2.0');
});

test('SKETCH-2827: a long symbol fades out on the last-used button', async ({ page }) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await chooseFromMonomerDatabase(page, 'LongPeptideSymbol', { search: 'long' });
  const button = await lastUsed(page);
  expect(button.text).toBe('LongPeptideSymbol');
  expect(button.toolTip).toBe('Long symbol glycine analog');
  await screenshotRegion(page, 'monomer-db-last-used-long-symbol.png', button);
});

test('SKETCH-2827: each page restores its own database monomer through the shared last-used button', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await chooseFromMonomerDatabase(page, 'Hyp');
  await clickWidget(page, 'nucleic_monomer_btn');
  // the other page has no last-used database monomer yet
  expect((await lastUsed(page)).text).toBe('');
  expect((await lastUsed(page)).enabled).toBe(false);
  await chooseFromMonomerDatabase(page, 'meA');
  for (const [pageButton, symbol] of [
    ['amino_monomer_btn', 'Hyp'],
    ['nucleic_monomer_btn', 'meA'],
    ['amino_monomer_btn', 'Hyp'],
  ]) {
    await clickWidget(page, pageButton);
    await expect.poll(async () => (await lastUsed(page)).text).toBe(symbol);
    expect((await lastUsed(page)).checked).toBe(true);
  }
  await place(page, -100);
  await clickWidget(page, 'nucleic_monomer_btn');
  await place(page, 100);
  expect((await waitAtoms(page, 2)).map((atom) => atom.label)).toEqual(['Hyp', 'meA']);
});

test('SKETCH-2827: the search is reset when the popup opens on a different page than it closed on', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabase(page);
  await searchMonomerDatabase(page, 'CHM');
  await closeMonomerDatabase(page);
  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('CHM');
  await closeMonomerDatabase(page);

  await clickWidget(page, 'nucleic_monomer_btn');
  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('');
  // the nucleic acid page's default filter shows nucleic acids and structures
  expect(await monomerDatabaseSymbols(page)).toEqual(['CHM1', 'CHM2', 'meA', 'NAx']);
  await searchMonomerDatabase(page, 'me');
  await closeMonomerDatabase(page);
  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('me');
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: typing in the search box does not trigger shortcuts, which work again once the popup closes', async ({
  page,
}) => {
  // Ala is already the active tool. (Clicking its checked button would open
  // its analog popup instead.)
  expect((await widgetState(page, 'ala_btn')).checked).toBe(true);
  await openMonomerDatabase(page);
  await page.keyboard.type('g', { delay: 10 });
  expect(await searchText(page)).toBe('g');
  expect((await widgetState(page, 'ala_btn')).checked).toBe(true);
  expect((await widgetState(page, 'gly_btn')).checked).toBe(false);
  await closeMonomerDatabase(page);
  await page.keyboard.press('g');
  await expect
    .poll(async () => (await widgetState(page, 'gly_btn')).checked, { timeout: 5000 })
    .toBe(true);
});

test('SKETCH-2827: Ctrl+A in the search box selects the search text and not the canvas contents', async ({
  page,
}) => {
  await pasteText(page, 'PEPTIDE1{A.G}$$$$V2.0');
  await openMonomerDatabase(page);
  await page.keyboard.type('xyz', { delay: 10 });
  await page.keyboard.press('ControlOrMeta+a');
  // typing replaces the selected search text
  await page.keyboard.type('am', { delay: 10 });
  expect(await searchText(page)).toBe('am');
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  await closeMonomerDatabase(page);
  // The fit button is restyled whenever there is a selection
  expect((await widgetState(page, 'fit_btn')).styleSheet).toBe('');
});

/**
 * Shrink the browser viewport and wait for Qt to see the new canvas size
 */
async function setViewportSize(page, width, height) {
  const before = await requireRect(page, 'widget:view');
  await page.setViewportSize({ width, height });
  await expect
    .poll(async () => {
      const view = await getRect(page, 'widget:view');
      return view !== null && (view.width !== before.width || view.height !== before.height);
    })
    .toBe(true);
}

function contains(outer, inner) {
  return (
    inner.x >= outer.x &&
    inner.y >= outer.y &&
    inner.x + inner.width <= outer.x + outer.width &&
    inner.y + inner.height <= outer.y + outer.height
  );
}

test('SKETCH-2920: the popup shrinks its table to stay inside a small canvas, including while it is open', async ({
  page,
}) => {
  await loadDatabase(page, 'chem_monomers.json');
  await openMonomerDatabase(page);
  const full = await requireRect(page, 'widget:monomer_db_popup');
  await closeMonomerDatabase(page);

  // Leave less room below the button than the popup needs. If the side bar
  // also leaves too little room above it, the table has to shrink.
  const button = await requireRect(page, 'widget:monomer_db_btn');
  const width = page.viewportSize().width;
  const roomBelow = full.height - 100;
  await setViewportSize(page, width, button.y + button.height + roomBelow);
  const shrinks = button.y < full.height;
  const checkInside = async () => {
    const size = page.viewportSize();
    const canvas = { x: 0, y: 0, width: size.width, height: size.height };
    const popup = await requireRect(page, 'widget:monomer_db_popup');
    expect(contains(canvas, popup)).toBe(true);
    for (const name of [
      'monomer_db_title_lbl',
      'monomer_db_close_btn',
      'monomer_db_search_le',
      'monomer_db_close_link',
    ]) {
      expect(contains(popup, await requireRect(page, `widget:${name}`)), name).toBe(true);
    }
    return popup;
  };
  await openMonomerDatabase(page);
  const fitted = await checkInside();
  if (shrinks) {
    expect(fitted.height).toBeLessThan(full.height);
  } else {
    expect(fitted.height).toBe(full.height);
  }
  // the table still scrolls to every row
  const { rowCount } = await monomerDatabaseTable(page);
  const last = await monomerDatabaseRow(page, rowCount - 1);
  expect(contains(fitted, last)).toBe(true);

  // shrink the canvas while the popup is open
  await setViewportSize(page, width, button.y + button.height + roomBelow - 40);
  if (shrinks) {
    await expect
      .poll(async () => (await requireRect(page, 'widget:monomer_db_popup')).height)
      .toBeLessThan(fitted.height);
  }
  await checkInside();
  await closeMonomerDatabase(page);
});

test('SKETCH-2920: below its minimum layout size, the popup keeps its top-left corner inside the canvas', async ({
  page,
}) => {
  const button = await requireRect(page, 'widget:monomer_db_btn');
  await setViewportSize(page, 300, button.y + button.height + 60);
  await openMonomerDatabase(page);
  const popup = await requireRect(page, 'widget:monomer_db_popup');
  const size = page.viewportSize();
  expect(popup.x).toBeGreaterThanOrEqual(0);
  expect(popup.y).toBeGreaterThanOrEqual(0);
  // the unavoidable overflow is to the right and bottom
  expect(popup.x + popup.width > size.width || popup.y + popup.height > size.height).toBe(true);
  // The x and Close are on the right, so they may be outside the canvas
  await page.keyboard.press('Escape');
  await expect.poll(() => isMonomerDatabaseOpen(page)).toBe(false);
});

test('SKETCH-2827: with a canvas selection, browsing never mutates and each click, Enter and Space commits once', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await pasteText(page, 'PEPTIDE1{A}$$$$V2.0');
  await selectAll(page);
  await openMonomerDatabase(page);
  // Tab moves the focus from the search box to the table, which makes the
  // first row current without selecting or committing it
  await page.keyboard.press('Tab');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('ArrowDown');
  let table = await monomerDatabaseTable(page);
  expect(table.currentRow).toBeGreaterThanOrEqual(0);
  expect(table.selectedRow).toBe(-1);
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{A}$$$$V2.0');
  expect((await lastUsed(page)).text).toBe('');

  // a click commits the clicked row
  const aib = await findMonomerDatabaseRow(page, 'Aib');
  await clickSelector(page, `table_row:${MONOMER_DB_TABLE}:${aib}`);
  await expect.poll(() => getExportedHelm(page)).toBe('PEPTIDE1{[Aib]}$$$$V2.0');
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  table = await monomerDatabaseTable(page);
  expect(table.currentRow).toBe(aib);
  expect(table.selectedRow).toBe(-1);
  expect((await lastUsed(page)).text).toBe('Aib');
  await screenshotRegion(page, 'monomer-db-mutation-current-row.png', table);

  // Enter and Space commit the current row. The rows are ac, Aib, am, Bz, ...
  for (const [move, key, symbol] of [
    ['ArrowDown', 'Enter', 'am'],
    ['ArrowDown', 'Space', 'Bz'],
  ]) {
    await page.keyboard.press(move);
    expect(await getExportedHelm(page)).not.toContain(`[${symbol}]`);
    await page.keyboard.press(key);
    await expect.poll(() => getExportedHelm(page)).toBe(`PEPTIDE1{[${symbol}]}$$$$V2.0`);
    expect(await isMonomerDatabaseOpen(page)).toBe(true);
    table = await monomerDatabaseTable(page);
    expect(table.selectedRow).toBe(-1);
    expect((await monomerDatabaseRow(page, table.currentRow)).cells[0]).toBe(symbol);
  }
  await closeMonomerDatabase(page);

  // each commit added exactly one undo entry
  for (const helm of ['PEPTIDE1{[am]}$$$$V2.0', 'PEPTIDE1{[Aib]}$$$$V2.0', 'PEPTIDE1{A}$$$$V2.0']) {
    await undo(page);
    await expect.poll(() => getExportedHelm(page)).toBe(helm);
  }
});

test('SKETCH-2827: Enter in tool mode activates the current row without closing the popup', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await chooseFromMonomerDatabase(page, 'Aib');
  await clickWidget(page, 'ala_btn');
  // reopening with a standard tool active leaves no row current or selected
  await openMonomerDatabase(page);
  let table = await monomerDatabaseTable(page);
  expect(table.currentRow).toBe(-1);
  expect(table.selectedRow).toBe(-1);
  // The search text "Aib" is kept, so Aib is the only row. Tab moves the focus
  // from the search box to the table, which makes it current without
  // selecting or activating it.
  await page.keyboard.press('Tab');
  table = await monomerDatabaseTable(page);
  expect(table.selectedRow).toBe(-1);
  expect((await lastUsed(page)).checked).toBe(false);
  expect((await monomerDatabaseRow(page, table.currentRow)).cells[0]).toBe('Aib');
  await page.keyboard.press('Enter');
  await expect.poll(async () => (await lastUsed(page)).checked).toBe(true);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  table = await monomerDatabaseTable(page);
  expect((await monomerDatabaseRow(page, table.selectedRow)).cells[0]).toBe('Aib');
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: removing a page’s database monomer while it is inactive restores the page’s standard tool', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  // inactive page
  await clickWidget(page, 'nucleic_monomer_btn');
  await chooseFromMonomerDatabase(page, 'meA');
  await clickWidget(page, 'amino_monomer_btn');
  // chem_monomers.json has no nucleic acids
  await loadDatabase(page, 'chem_monomers.json');
  await clickWidget(page, 'nucleic_monomer_btn');
  await expect.poll(async () => (await widgetState(page, 'na_rna_btn')).checked).toBe(true);
  expect((await lastUsed(page)).text).toBe('');
  expect((await lastUsed(page)).checked).toBe(false);

  // atomistic mode
  await clickWidget(page, 'amino_monomer_btn');
  await chooseFromMonomerDatabase(page, 'CHM1');
  await clickWidget(page, 'atomistic_btn');
  // na_monomers.json has no structures
  await loadDatabase(page, 'na_monomers.json');
  await clickWidget(page, 'monomeric_btn');
  await expect.poll(async () => (await widgetState(page, 'ala_btn')).checked).toBe(true);
  expect((await lastUsed(page)).text).toBe('');
  await place(page);
  expect(await getExportedHelm(page)).toBe('PEPTIDE1{A}$$$$V2.0');
});

test('SKETCH-2827: a structure can be chosen from either page', async ({ page }) => {
  await loadDatabase(page, 'chem_monomers.json');
  await chooseFromMonomerDatabase(page, 'CHM2');
  await place(page);
  expect(await getExportedHelm(page)).toContain('CHEM1{[CHM2]}');
  await clear(page);
  await clickWidget(page, 'nucleic_monomer_btn');
  await chooseFromMonomerDatabase(page, 'CHM1');
  expect((await lastUsed(page)).text).toBe('CHM1');
  await place(page);
  expect(await getExportedHelm(page)).toContain('CHEM1{[CHM1]}');
  // each page remembers its own structure
  await clickWidget(page, 'amino_monomer_btn');
  await expect.poll(async () => (await lastUsed(page)).text).toBe('CHM2');
});
