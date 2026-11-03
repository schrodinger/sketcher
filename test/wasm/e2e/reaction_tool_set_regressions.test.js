import { expect, test } from '@playwright/test';
import {
  center,
  clickWidget,
  isSketcherEmpty,
  isWidgetVisible,
  mouseClick,
  openSketcher,
  redo,
  undo,
  waitAtoms,
  widgetState,
} from './e2e_helpers.js';

const OBJECT_POSITION = [0, 0];
const ATOM_POSITION = [150, 0];

/**
 * Wait for the workspace to report whether it is empty, and for the side bar's
 * tool set buttons to agree: an atomistic tool set is always available, but
 * the monomeric one only while the workspace is empty.
 */
async function expectWorkspace(page, empty) {
  await expect.poll(() => isSketcherEmpty(page)).toBe(empty);
  await expect.poll(async () => (await widgetState(page, 'monomeric_btn')).enabled).toBe(empty);
  expect((await widgetState(page, 'atomistic_btn')).enabled).toBe(true);
}

async function expectToolSet(page, monomeric) {
  await expect.poll(async () => (await widgetState(page, 'monomeric_btn')).checked).toBe(monomeric);
  expect((await widgetState(page, 'atomistic_btn')).checked).toBe(!monomeric);
  await expect.poll(() => isWidgetVisible(page, 'amino_monomer_btn')).toBe(monomeric);
}

// SKETCH-2925: reaction arrows and plus signs aren't allowed in monomeric
// mode, so drawing one must keep the user from switching to monomeric tools.
for (const [tool, description] of [
  ['rxn_arrow', 'reaction arrow'],
  ['rxn_plus', 'reaction plus sign'],
]) {
  test.describe(`SKETCH-2925 ${description}`, () => {
    test(`a ${description} prevents switching to monomeric tools`, async ({ page }) => {
      const sk = await openSketcher(page);
      await expectWorkspace(page, true);
      await expectToolSet(page, false);

      await sk.click_tool(tool);
      await sk.click_sketcher(...OBJECT_POSITION);
      await expectWorkspace(page, false);

      // Click where the user would. clickWidget waits for a control to become
      // enabled, so it can't click a disabled one.
      const monomeric_btn = center(await widgetState(page, 'monomeric_btn'));
      await mouseClick(page, monomeric_btn.x, monomeric_btn.y);
      await expectToolSet(page, false);

      // Undoing the drawing makes monomeric tools available again, and redoing
      // it while they're shown returns to atomistic tools.
      await undo(page);
      await expectWorkspace(page, true);
      await clickWidget(page, 'monomeric_btn');
      await expectToolSet(page, true);
      await redo(page);
      await expectWorkspace(page, false);
      await expectToolSet(page, false);
    });

    test(`erasing atoms but not a ${description} keeps monomeric tools unavailable`, async ({
      page,
    }) => {
      const sk = await openSketcher(page);
      await sk.click_tool(tool);
      await sk.click_sketcher(...OBJECT_POSITION);
      await sk.click_tool('C');
      await sk.click_sketcher(...ATOM_POSITION);
      await waitAtoms(page, 1);
      await expectWorkspace(page, false);

      await sk.click_tool('erase');
      await sk.click_sketcher(...ATOM_POSITION);
      await waitAtoms(page, 0);
      await expectWorkspace(page, false);

      await sk.click_sketcher(...OBJECT_POSITION);
      await expectWorkspace(page, true);
      await undo(page);
      await expectWorkspace(page, false);
      await redo(page);
      await expectWorkspace(page, true);

      await undo(page);
      await expectWorkspace(page, false);
      await sk.click_button('clear');
      await expectWorkspace(page, true);
      await clickWidget(page, 'monomeric_btn');
      await expectToolSet(page, true);
    });
  });
}
