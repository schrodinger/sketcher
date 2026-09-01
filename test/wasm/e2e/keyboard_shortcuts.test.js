import { expect, test } from '@playwright/test';
import {
  waitForSketcherReady,
  focusCanvas,
  getCanvasCenter,
  getExportedSmiles,
} from './e2e_helpers.js';

// Ported (behavioral subset) from `suite_2D_sketcher_new/tst_hidden_shortcuts`.
// The original asserted that typing an element key checks the matching draw-tool
// button (button state) and then renders; it covered f,h,n,o,p,s,i,b,k,u plus
// bond-order digits. Querying Qt button `checked`/`enabled` state needs a WASM
// binding the standalone harness doesn't expose yet, so here we verify the
// user-observable half: typing an element key selects that draw tool, so a
// subsequent click on the empty canvas places that atom.
//
// Limited to elements whose lone-atom SMILES is unambiguous. n/o/s are already
// covered in drawing.test.js; the bond-order digit shortcuts and tool-state
// assertions are deferred until atom/bond hit-testing + button-state bindings
// exist.

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Keyboard element shortcuts', () => {
  for (const { key, smiles } of [
    { key: 'f', smiles: 'F' },
    { key: 'p', smiles: 'P' },
    { key: 'i', smiles: 'I' },
  ]) {
    test(`pressing "${key}" then clicking places ${smiles}`, async ({ page }) => {
      const center = await getCanvasCenter(page);

      await focusCanvas(page);
      await page.keyboard.press(key);
      await page.mouse.click(center.x, center.y);

      await expect.poll(() => getExportedSmiles(page), { timeout: 5000 }).toBe(smiles);
    });
  }
});
