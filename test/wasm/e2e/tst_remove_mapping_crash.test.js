import { expect, test } from '@playwright/test';
import {
  canvasPoint,
  clickItem,
  clickPopupTool,
  clickWidget,
  getExportedSmiles,
  loadStructure,
  mouseClick,
  openSketcher,
  waitAtoms,
} from './e2e_helpers.js';

const REACTION = '[CH3:1][CH3:2]>>[CH3:1][CH3:2]';
// Each iteration gives the freed atom item and atoms another chance to be
// reused, which is what turns the stale read into a crash or a wrong edit
const ITERATIONS = 20;

test.setTimeout(180_000);

/**
 * The remove mapping tool used to keep a pointer to the atom item from the
 * previous click.  Removing a mapping deletes all atom items, and undo then
 * replaces the molecule, so a later click on empty space read both freed
 * objects.  Depending on how that memory had been reused, this either crashed
 * or removed the wrong mapping and broke later edits.
 */
test.describe('tst_remove_mapping_crash', () => {
  test('click empty space after removing a mapping and undoing', async ({ page }) => {
    let browserCrashed = false;
    page.on('crash', () => {
      browserCrashed = true;
    });
    await openSketcher(page);
    await loadStructure(page, REACTION);
    await waitAtoms(page, 4);
    await clickPopupTool(page, 'remove_mapping_btn');
    // Below the horizontal reaction, away from all atoms and the arrow
    const empty = await canvasPoint(page, 0, 150);

    for (let i = 0; i < ITERATIONS; i += 1) {
      const message = `iteration ${i}`;
      // atom 0 is the reactant mapped to 1, so this removes mapping 1 from
      // the reactant and the product
      await clickItem(page, 'atom', 0);
      await expect.poll(() => getExportedSmiles(page), { message }).not.toContain(':1]');

      await clickWidget(page, 'undo_btn');
      await expect.poll(() => getExportedSmiles(page), { message }).toContain(':1]');
      // exporting also allocates, which makes reuse of the freed memory likely
      const smiles = await getExportedSmiles(page);

      await mouseClick(page, empty.x, empty.y);
      expect(browserCrashed, message).toBe(false);
      expect(await getExportedSmiles(page), message).toBe(smiles);
    }
  });
});
