import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getCanvasCenter,
  getDrawingAreaCenter,
  getExportedSmiles,
  isSketcherEmpty,
  waitForMoleculeChange,
  waitForRender,
  modifierKey,
  clickToolbarButton,
} from './e2e_helpers.js';

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('E2E Interaction Tests', () => {
  test('click on empty canvas places a carbon atom', async ({ page }) => {
    const center = await getCanvasCenter(page);

    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));

    expect(await isSketcherEmpty(page)).toBe(false);
    expect(await getExportedSmiles(page)).toBe('C');
  });

  // Each element shortcut gets its own test with a fresh page to avoid
  // interference from programmatic clears, which can leave the sketcher
  // in a transient state that ignores subsequent mouse clicks.
  for (const { key, smiles: expected } of [
    { key: 'N', smiles: 'N' },
    { key: 'O', smiles: 'O' },
    { key: 'S', smiles: 'S' },
  ]) {
    test(`keyboard shortcut "${key}" places ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await focusCanvas(page);
      await page.keyboard.press(key);
      await page.mouse.click(center.x, center.y);
      await page.waitForFunction(() => !Module.sketcher_is_empty(), { timeout: 5000 });

      expect(await getExportedSmiles(page)).toBe(expected);
    });
  }

  test('click-drag draws a bond', async ({ page }) => {
    const center = await getCanvasCenter(page);
    const dragDistance = 100;

    await waitForMoleculeChange(page, async () => {
      await page.mouse.move(center.x, center.y);
      await page.mouse.down();
      await page.mouse.move(center.x + dragDistance, center.y, { steps: 10 });
      await page.mouse.up();
    });

    expect(await isSketcherEmpty(page)).toBe(false);
    expect(await getExportedSmiles(page)).toBe('CC');
  });

  test('undo and redo', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place an atom
    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));
    expect(await isSketcherEmpty(page)).toBe(false);
    expect(await getExportedSmiles(page)).toBe('C');

    // Undo — sketcher should be empty again
    const mod = modifierKey();
    await focusCanvas(page);
    await waitForMoleculeChange(page, () => page.keyboard.press(`${mod}+z`));
    expect(await isSketcherEmpty(page)).toBe(true);

    // Redo — atom should reappear
    await focusCanvas(page);
    await waitForMoleculeChange(page, () => page.keyboard.press(`${mod}+Shift+z`));
    expect(await isSketcherEmpty(page)).toBe(false);
    expect(await getExportedSmiles(page)).toBe('C');
  });

  // Toolbar button click tests — each gets a fresh page
  for (const { button, smiles: expected } of [
    { button: 'c_btn', smiles: 'C' },
    { button: 'n_btn', smiles: 'N' },
    { button: 'o_btn', smiles: 'O' },
    { button: 's_btn', smiles: 'S' },
    { button: 'f_btn', smiles: 'F' },
    { button: 'p_btn', smiles: 'P' },
  ]) {
    test(`toolbar button "${button}" places ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      // Click the element button on the toolbar to select the tool
      await clickToolbarButton(page, button);
      // Click on the canvas to place the atom
      await page.mouse.click(center.x, center.y);
      await page.waitForFunction(() => !Module.sketcher_is_empty(), { timeout: 5000 });

      expect(await getExportedSmiles(page)).toBe(expected);
    });
  }

  // Bond tool tests — select tool, then click-drag to draw
  for (const { button, smiles: expected } of [{ button: 'bond_order_btn', smiles: 'C=C' }]) {
    test(`bond tool "${button}" draws ${expected}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await clickToolbarButton(page, button);
      await waitForMoleculeChange(page, async () => {
        await page.mouse.move(center.x, center.y);
        await page.mouse.down();
        await page.mouse.move(center.x + 100, center.y, { steps: 10 });
        await page.mouse.up();
      });

      expect(await getExportedSmiles(page)).toBe(expected);
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

      await clickToolbarButton(page, button);
      await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));

      expect(await getExportedSmiles(page)).toBe(expected);
    });
  }

  test('erase tool removes an atom', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place an atom
    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));
    expect(await isSketcherEmpty(page)).toBe(false);

    // Switch to erase tool and click the atom
    await clickToolbarButton(page, 'erase_btn');
    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));
    expect(await isSketcherEmpty(page)).toBe(true);
  });

  test('increase charge tool adds positive charge', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Place a carbon atom
    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));
    expect(await getExportedSmiles(page)).toBe('C');

    // Switch to increase charge tool and click the atom
    await clickToolbarButton(page, 'increase_charge_btn');
    await waitForMoleculeChange(page, () => page.mouse.click(center.x, center.y));

    const smiles = await getExportedSmiles(page);
    expect(smiles).toContain('+');
  });

  test('import molecule via API then extend by drawing', async ({ page }) => {
    // Import a molecule programmatically (the LiveDesign workflow)
    await page.evaluate(() => {
      Module.sketcher_import_text('C');
    });
    expect(await getExportedSmiles(page)).toBe('C');

    // Fit the view so the imported atom is centered in the drawing area
    await clickToolbarButton(page, 'fit_btn');
    await waitForRender(page);

    // Use the drawing area center (not canvas center) to hit the atom
    const center = await getDrawingAreaCenter(page);

    // Draw a bond from the atom outward to extend the molecule
    await clickToolbarButton(page, 'single_bond_btn');
    await waitForMoleculeChange(page, async () => {
      await page.mouse.move(center.x, center.y);
      await page.mouse.down();
      await page.mouse.move(center.x + 100, center.y, { steps: 10 });
      await page.mouse.up();
    });

    const smiles = await getExportedSmiles(page);
    expect(smiles).toBe('CC');
  });

  test('visual regression snapshot', async ({ page }) => {
    const center = await getCanvasCenter(page);

    // Draw a two-atom molecule via click-drag
    await waitForMoleculeChange(page, async () => {
      await page.mouse.move(center.x, center.y);
      await page.mouse.down();
      await page.mouse.move(center.x + 100, center.y, { steps: 10 });
      await page.mouse.up();
    });

    // Move the mouse away so the cursor doesn't affect the screenshot
    await page.mouse.move(0, 0);
    await waitForRender(page, 300);

    const canvas = page.locator('#screen canvas');
    expect(await canvas.screenshot()).toMatchSnapshot('e2e-draw-bond.png');
  });
});
