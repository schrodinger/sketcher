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
  chooseFilterItem,
  chooseFromMonomerDatabase,
  closeMonomerDatabase,
  closeMonomerDatabaseFilter,
  closePolymerTypeDropdown,
  findMonomerDatabaseRow,
  isMonomerDatabaseFilterOpen,
  isMonomerDatabaseOpen,
  loadDatabase,
  monomerDatabaseRow,
  monomerDatabaseSymbols,
  monomerDatabaseTable,
  openComboDropdown,
  openMonomerDatabase,
  openMonomerDatabaseFilter,
  openMonomerSketcher,
  searchMonomerDatabase,
  setPolymerTypeChecked,
  tabFromSearchToMonomerDatabaseTable,
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

/** The Monomer Database table's names, lowercased, in display order */
async function monomerDatabaseNames(page) {
  const { rowCount } = await monomerDatabaseTable(page);
  const names = [];
  for (let row = 0; row < rowCount; row += 1) {
    names.push((await monomerDatabaseRow(page, row)).cells[1].toLowerCase());
  }
  return names;
}

/**
 * Click the Monomer Database table's Name header twice, which sorts by name,
 * descending, when the table is sorted by another column or by name,
 * descending.
 */
async function sortMonomerDatabaseByNameDescending(page) {
  const table = await monomerDatabaseTable(page);
  const x = table.x + table.width / 2;
  const y = table.y + 10;
  await mouseClick(page, x, y);
  // Click somewhere else in the header, so that Qt doesn't treat the second
  // click as a double-click, which wouldn't change the sort again
  await mouseClick(page, x + 20, y);
  const names = await monomerDatabaseNames(page);
  expect(names).toEqual([...names].sort().reverse());
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
  await sortMonomerDatabaseByNameDescending(page);
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
  // Moving the focus from the search box to the table makes the first row
  // current without selecting or committing it
  await tabFromSearchToMonomerDatabaseTable(page);
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
  // The search text "Aib" is kept, so Aib is the only row. Moving the focus
  // from the search box to the table makes it current without selecting or
  // activating it.
  await tabFromSearchToMonomerDatabaseTable(page);
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

/** The filter button's style sheet shows the orange icon when filtering */
async function isFilterButtonFiltering(page) {
  return (await widgetState(page, 'monomer_db_filter_btn')).styleSheet.includes(
    'filter-menu_lb_filtering.png',
  );
}

/** Return a filter popup combo box's displayed text, enabled state and boldness */
async function comboState(page, combo) {
  const { text, enabled, bold } = await widgetState(page, combo);
  return { text, enabled, bold };
}

async function filtersApplied(page) {
  return (await widgetState(page, 'filters_applied_lbl')).text;
}

/** Move the mouse off the filter button so that it isn't drawn hovered */
async function moveMouseToPopupTitle(page) {
  const title = await requireRect(page, 'widget:monomer_db_title_lbl');
  await page.mouse.move(title.x + 5, title.y + title.height / 2);
}

test('SKETCH-2827: polymer types are toggled in a dropdown nested in two popups, and filter the table live', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabaseFilter(page);
  // the amino acid page's default
  expect(await comboState(page, 'polymer_type_combo')).toEqual({
    text: 'Peptide, Structure',
    enabled: true,
    bold: true,
  });
  expect((await comboState(page, 'monomer_type_combo')).enabled).toBe(true);
  expect((await comboState(page, 'natural_analog_combo')).enabled).toBe(true);
  expect(await filtersApplied(page)).toBe('1 filter applied');
  expect((await monomerDatabaseTable(page)).headers).toEqual(['Symbol', 'Name', 'Type', 'Analog']);

  // checking Nucleic acid keeps the dropdown, the filter popup and the
  // Monomer Database popup open
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  expect((await widgetState(page, 'polymer_type_combo')).open).toBe(true);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(true);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  expect(await comboState(page, 'polymer_type_combo')).toEqual({
    text: 'Any',
    enabled: true,
    bold: false,
  });
  expect(await comboState(page, 'monomer_type_combo')).toEqual({
    text: 'Any',
    enabled: false,
    bold: false,
  });
  expect(await comboState(page, 'natural_analog_combo')).toEqual({
    text: 'Any or None',
    enabled: false,
    bold: false,
  });
  expect(await filtersApplied(page)).toBe('No filters applied');
  expect(await matchCount(page)).toBe('11 matches');

  // nucleic acids only
  await setPolymerTypeChecked(page, 'Peptide', false);
  await setPolymerTypeChecked(page, 'Structure', false);
  expect((await comboState(page, 'polymer_type_combo')).text).toBe('Nucleic acid');
  expect((await comboState(page, 'monomer_type_combo')).enabled).toBe(false);
  expect((await comboState(page, 'natural_analog_combo')).enabled).toBe(true);
  expect(await monomerDatabaseSymbols(page)).toEqual(['meA', 'NAx']);
  expect(await matchCount(page)).toBe('2 matches');
  // a single polymer type hides the Type column
  expect((await monomerDatabaseTable(page)).headers).toEqual(['Symbol', 'Name', 'Analog']);

  // unchecking the last polymer type empties the table and keeps the dropdown
  // open, and checking it again restores the rows
  await setPolymerTypeChecked(page, 'Nucleic acid', false);
  expect(await comboState(page, 'polymer_type_combo')).toEqual({
    text: 'None',
    enabled: true,
    bold: true,
  });
  expect((await widgetState(page, 'polymer_type_combo')).open).toBe(true);
  expect((await comboState(page, 'natural_analog_combo')).enabled).toBe(false);
  expect(await monomerDatabaseSymbols(page)).toEqual([]);
  expect(await matchCount(page)).toBe('0 matches');
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  expect(await monomerDatabaseSymbols(page)).toEqual(['meA', 'NAx']);

  // Escape closes only the dropdown
  await closePolymerTypeDropdown(page);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(true);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);

  // the analog filter applies to nucleic acids
  await chooseFilterItem(page, 'natural_analog_combo', 'A');
  expect((await comboState(page, 'natural_analog_combo')).bold).toBe(true);
  expect(await filtersApplied(page)).toBe('2 filters applied');
  expect(await monomerDatabaseSymbols(page)).toEqual(['meA']);
  await closeMonomerDatabaseFilter(page);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: after an item is clicked, a click elsewhere in the filter popup closes only the polymer type dropdown', async ({
  page,
}) => {
  await openMonomerDatabaseFilter(page);
  // clicking an item must not leave the filter popup unable to close the
  // dropdown
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  await clickWidget(page, 'polymer_type_lbl');
  await expect.poll(async () => (await widgetState(page, 'polymer_type_combo')).open).toBe(false);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(true);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);

  // a reopened dropdown is still the popup that Escape closes
  await openComboDropdown(page, 'polymer_type_combo');
  await page.keyboard.press('Escape');
  await expect.poll(async () => (await widgetState(page, 'polymer_type_combo')).open).toBe(false);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(true);
  await closeMonomerDatabaseFilter(page);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: after a filter combo box is used, a click on the Monomer Database popup still closes the filter popup', async ({
  page,
}) => {
  /** Click the Monomer Database popup outside of the filter popup */
  async function clickMonomerDatabasePopup() {
    await clickWidget(page, 'monomer_db_title_lbl');
  }

  // a dropdown closed by a click in the filter popup
  await openMonomerDatabaseFilter(page);
  await openComboDropdown(page, 'monomer_type_combo');
  await clickWidget(page, 'polymer_type_lbl');
  await expect.poll(async () => (await widgetState(page, 'monomer_type_combo')).open).toBe(false);
  await clickMonomerDatabasePopup();
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);

  // a dropdown closed by choosing an item
  await openMonomerDatabaseFilter(page);
  await chooseFilterItem(page, 'monomer_type_combo', 'Backbone');
  await clickMonomerDatabasePopup();
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);

  // an open dropdown takes one click to close, and the filter popup another
  await openMonomerDatabaseFilter(page);
  await openComboDropdown(page, 'polymer_type_combo');
  await clickMonomerDatabasePopup();
  await expect.poll(async () => (await widgetState(page, 'polymer_type_combo')).open).toBe(false);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(true);
  await clickMonomerDatabasePopup();
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: Backbone and an analog survive checking and unchecking Nucleic acid, but not switching to it', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabaseFilter(page);
  await setPolymerTypeChecked(page, 'Structure', false);
  await closePolymerTypeDropdown(page);
  await chooseFilterItem(page, 'monomer_type_combo', 'Backbone');
  await chooseFilterItem(page, 'natural_analog_combo', 'A');
  expect(await monomerDatabaseSymbols(page)).toEqual(['dAla']);
  expect(await filtersApplied(page)).toBe('3 filters applied');

  // checking Nucleic acid disables both combos, which show Any and Any or None
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  expect(await comboState(page, 'monomer_type_combo')).toEqual({
    text: 'Any',
    enabled: false,
    bold: false,
  });
  expect(await comboState(page, 'natural_analog_combo')).toEqual({
    text: 'Any or None',
    enabled: false,
    bold: false,
  });
  expect(await filtersApplied(page)).toBe('1 filter applied');
  expect(await matchCount(page)).toBe('9 matches');

  // unchecking it restores them
  await setPolymerTypeChecked(page, 'Nucleic acid', false);
  expect(await comboState(page, 'monomer_type_combo')).toEqual({
    text: 'Backbone',
    enabled: true,
    bold: true,
  });
  expect(await comboState(page, 'natural_analog_combo')).toEqual({
    text: 'A',
    enabled: true,
    bold: true,
  });
  expect(await monomerDatabaseSymbols(page)).toEqual(['dAla']);

  // checking Nucleic acid and then unchecking Peptide forgets both
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  await setPolymerTypeChecked(page, 'Peptide', false);
  expect((await comboState(page, 'natural_analog_combo')).text).toBe('Any or None');
  expect((await comboState(page, 'natural_analog_combo')).enabled).toBe(true);
  await setPolymerTypeChecked(page, 'Peptide', true);
  await setPolymerTypeChecked(page, 'Nucleic acid', false);
  expect(await comboState(page, 'monomer_type_combo')).toEqual({
    text: 'Any',
    enabled: true,
    bold: false,
  });
  expect(await comboState(page, 'natural_analog_combo')).toEqual({
    text: 'Any or None',
    enabled: true,
    bold: false,
  });
  await closePolymerTypeDropdown(page);
  await closeMonomerDatabaseFilter(page);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: Clear All shows every monomer on either page and keeps the search and sort', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  for (const pageButton of ['amino_monomer_btn', 'nucleic_monomer_btn']) {
    await clickWidget(page, pageButton);
    await openMonomerDatabase(page);
    await searchMonomerDatabase(page, 'a');
    await sortMonomerDatabaseByNameDescending(page);
    expect(await isFilterButtonFiltering(page)).toBe(true);
    if (pageButton === 'amino_monomer_btn') {
      await moveMouseToPopupTitle(page);
      await screenshotRegion(
        page,
        'monomer-db-filter-button-filtering.png',
        await requireRect(page, 'widget:monomer_db_filter_btn'),
      );
    }

    await openMonomerDatabaseFilter(page);
    await clickWidget(page, 'filter_clear_all_link');
    await expect.poll(() => filtersApplied(page)).toBe('No filters applied');
    expect(await comboState(page, 'polymer_type_combo')).toEqual({
      text: 'Any',
      enabled: true,
      bold: false,
    });
    expect((await comboState(page, 'monomer_type_combo')).text).toBe('Any');
    expect((await comboState(page, 'natural_analog_combo')).text).toBe('Any or None');
    expect(await isFilterButtonFiltering(page)).toBe(false);
    await closeMonomerDatabaseFilter(page);
    expect(await searchText(page)).toBe('a');
    // every polymer type, still sorted by name, descending
    const symbols = await monomerDatabaseSymbols(page);
    expect(symbols).toContain('meA');
    expect(symbols).toContain('dAla');
    const names = await monomerDatabaseNames(page);
    expect(names).toEqual([...names].sort().reverse());
    expect((await monomerDatabaseTable(page)).headers).toEqual([
      'Symbol',
      'Name',
      'Type',
      'Analog',
    ]);
    if (pageButton === 'amino_monomer_btn') {
      await moveMouseToPopupTitle(page);
      await screenshotRegion(
        page,
        'monomer-db-filter-button-normal.png',
        await requireRect(page, 'widget:monomer_db_filter_btn'),
      );
    }
    await closeMonomerDatabase(page);
  }
});

test('SKETCH-2827: clicking the filter button while the filter popup is open only closes it', async ({
  page,
}) => {
  await openMonomerDatabaseFilter(page);
  await clickWidget(page, 'monomer_db_filter_btn');
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  // give a reopened popup time to appear
  await page.waitForTimeout(500);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: Tab moves the focus from the search box to the filter button and then to the table', async ({
  page,
}) => {
  await openMonomerDatabase(page);
  // the popup opens with the focus in the search box, so one Tab reaches the
  // filter button, which Space presses
  await page.keyboard.press('Tab');
  await page.keyboard.press('Space');
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(true);
  // closing the filter popup returns the focus to the filter button
  await page.keyboard.press('Escape');
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  expect((await monomerDatabaseTable(page)).currentRow).toBe(-1);
  // the next Tab reaches the table, which makes its first row current
  await page.keyboard.press('Tab');
  await expect.poll(async () => (await monomerDatabaseTable(page)).currentRow).toBe(0);
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(false);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: updating the database keeps the filter, search and sort on the same page', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabaseFilter(page);
  await chooseFilterItem(page, 'monomer_type_combo', 'Terminal');
  expect(await monomerDatabaseSymbols(page)).toEqual(['ac', 'am', 'Bz']);
  // with the filter popup open, a click on the load link only closes the
  // filter popup
  await clickWidget(page, 'monomer_db_load_link');
  await expect.poll(() => isMonomerDatabaseFilterOpen(page)).toBe(false);
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  await searchMonomerDatabase(page, 'a');

  await setFilePickerResult(page, path.join(FIXTURES, 'chem_monomers.json'));
  await clickLoadLink(page, 'Update Database...');
  expect(await isMonomerDatabaseFilterOpen(page)).toBe(false);
  await expectDialog(page, 'Monomer Database', 30000);
  await closeDialog(page);
  expect(await isMonomerDatabaseOpen(page)).toBe(false);

  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('a');
  // chem_monomers.json has no terminal peptides
  expect(await monomerDatabaseSymbols(page)).toEqual(['ac', 'am']);
  await openMonomerDatabaseFilter(page);
  expect((await comboState(page, 'monomer_type_combo')).text).toBe('Terminal');
  await closeMonomerDatabaseFilter(page);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: choosing an RNA row on the amino acid page switches pages, and the filter is kept there', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await openMonomerDatabaseFilter(page);
  await setPolymerTypeChecked(page, 'Nucleic acid', true);
  await setPolymerTypeChecked(page, 'Structure', false);
  await closePolymerTypeDropdown(page);
  await closeMonomerDatabaseFilter(page);
  const before = await requireRect(page, 'widget:monomer_db_popup');

  await chooseFromMonomerDatabase(page, 'meA', { close: false });
  await expect
    .poll(async () => (await widgetState(page, 'nucleic_monomer_btn')).checked)
    .toBe(true);
  expect((await lastUsed(page)).text).toBe('meA');
  // the popup stays open where it was
  expect(await isMonomerDatabaseOpen(page)).toBe(true);
  const after = await requireRect(page, 'widget:monomer_db_popup');
  expect({ x: after.x, y: after.y }).toEqual({ x: before.x, y: before.y });
  await closeMonomerDatabase(page);

  await openMonomerDatabase(page);
  expect(await searchText(page)).toBe('meA');
  await openMonomerDatabaseFilter(page);
  expect((await comboState(page, 'polymer_type_combo')).text).toBe('Peptide, Nucleic acid');
  await closeMonomerDatabaseFilter(page);
  await closeMonomerDatabase(page);
});

test('SKETCH-2827: with a canvas selection, browsing to an RNA row keeps the page and committing it switches pages', async ({
  page,
}) => {
  await loadDatabase(page, 'mixed_monomers.json');
  await pasteText(page, 'PEPTIDE1{A}$$$$V2.0');
  await selectAll(page);
  for (const toolButton of ['select_tool_btn', 'move_rotate_btn']) {
    if (toolButton === 'move_rotate_btn') await clickWidget(page, toolButton);
    await expect.poll(async () => (await widgetState(page, toolButton)).checked).toBe(true);
    expect((await widgetState(page, 'amino_monomer_btn')).checked).toBe(true);

    await openMonomerDatabaseFilter(page);
    await setPolymerTypeChecked(page, 'Nucleic acid', true);
    await closePolymerTypeDropdown(page);
    await closeMonomerDatabaseFilter(page);
    await searchMonomerDatabase(page, 'meA');
    // moving the focus to the table makes meA current
    await tabFromSearchToMonomerDatabaseTable(page);
    const table = await monomerDatabaseTable(page);
    expect((await monomerDatabaseRow(page, table.currentRow)).cells[0]).toBe('meA');
    expect((await widgetState(page, 'amino_monomer_btn')).checked).toBe(true);

    await page.keyboard.press('Enter');
    await expect
      .poll(async () => (await widgetState(page, 'nucleic_monomer_btn')).checked)
      .toBe(true);
    expect((await widgetState(page, toolButton)).checked).toBe(true);
    // the selected peptide can't become a nucleic acid monomer
    expect(await getExportedHelm(page)).toBe('PEPTIDE1{A}$$$$V2.0');
    expect(await isMonomerDatabaseOpen(page)).toBe(true);
    await closeMonomerDatabase(page);
    await clickWidget(page, 'amino_monomer_btn');
  }
});
