import { expect, test } from '@playwright/test';
import {
  clickWidget,
  drawBond,
  drawElement,
  exportText,
  openSketcher,
} from '../wrappers/sketcher_wasm.js';

test.beforeEach(async ({ page }) => {
  await openSketcher(page);
});

test.describe('ported Squish tools and shortcuts', () => {
  // Source coverage: tst_tools and tst_hidden_shortcuts.
  for (const [key, expectedSmiles] of [
    ['N', 'N'],
    ['O', 'O'],
    ['S', 'S'],
  ]) {
    test(`keyboard tool ${key} places ${expectedSmiles}`, async ({ page }) => {
      await drawElement(page, key);
      await expect.poll(() => exportText(page)).toBe(expectedSmiles);
    });
  }

  test('element toolbar button selects nitrogen', async ({ page }) => {
    await clickWidget(page, 'n_btn');
    await drawBond(page);
    await expect.poll(() => exportText(page)).toBe('NN');
  });

  test('double-bond toolbar tool draws a double bond', async ({ page }) => {
    await clickWidget(page, 'bond_order_btn');
    await drawBond(page);
    await expect.poll(() => exportText(page)).toBe('C=C');
  });
});
