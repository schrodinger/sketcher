import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  getCanvasCenter,
  getExportedSmiles,
  isSketcherEmpty,
  clickWidget,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Toolbar', () => {
  // Element button tests — each gets a fresh page
  for (const { button, smiles: expected } of [
    { button: 'c_btn', smiles: 'C' },
    { button: 'n_btn', smiles: 'N' },
    { button: 'o_btn', smiles: 'O' },
    { button: 's_btn', smiles: 'S' },
    { button: 'f_btn', smiles: 'F' },
    { button: 'p_btn', smiles: 'P' },
  ]) {
    test(`element button "${button}" places ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await clickWidget(page, button);
      await page.mouse.click(center.x, center.y);

      await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe(expected);
    });
  }

  // Bond tool tests — select tool, then click-drag to draw
  for (const { button, smiles: expected } of [{ button: 'bond_order_btn', smiles: 'C=C' }]) {
    test(`bond tool "${button}" draws ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await clickWidget(page, button);
      await page.mouse.move(center.x, center.y);
      await page.mouse.down();
      await page.mouse.move(center.x + 100, center.y, { steps: 10 });
      await page.mouse.up();

      await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe(expected);
    });
  }

  // Ring tool tests — select tool, then click to place
  for (const { button, smiles: expected } of [
    { button: 'cyclohexane_btn', smiles: 'C1CCCCC1' },
    { button: 'benzene_btn', smiles: 'C1=CC=CC=C1' },
    { button: 'cyclopentane_btn', smiles: 'C1CCCC1' },
  ]) {
    test(`ring tool "${button}" places ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await clickWidget(page, button);
      await page.mouse.click(center.x, center.y);

      await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe(expected);
    });
  }

  test('erase tool removes an atom', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place an atom
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => isSketcherEmpty(page), { timeout: 5000 }).toBe(false);

    // Switch to erase tool and click the atom
    await clickWidget(page, 'erase_btn');
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => isSketcherEmpty(page), { timeout: 5000 }).toBe(true);
  });

  test('increase charge tool adds positive charge', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place a carbon atom
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe('C');

    // Switch to increase charge tool and click the atom
    await clickWidget(page, 'increase_charge_btn');
    await page.mouse.click(center.x, center.y);
    await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toContain('+');
  });
});
