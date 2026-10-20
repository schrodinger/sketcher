import { expect, test } from '@playwright/test';
import { fileURLToPath } from 'node:url';
import { clickPopupTool, holdPopup, waitAtoms } from './e2e_helpers.js';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import {
  clickWidget,
  drawingAreaCenter,
  exportText,
  mouseClick,
  mouseDrag,
  widgetState,
  widgetRect,
} from './sketcher/wrappers/sketcher_wasm.js';

// SKETCH-2911: use the case's ChEMBL database, including a long symbol with
// both backbone attachment sites so the drag can create a connected pair.
const MONOMER_DATABASE = fileURLToPath(
  new URL('./fixtures/chembl_monomers_no_id.json', import.meta.url),
);
const LONG_SYMBOL = 'dSer(tBu)';
const ORIGINAL_MONOMERS = 6;

async function peptideSequences(page) {
  const helm = await exportText(page, 'HELM');
  return [...helm.split('$')[0].matchAll(/PEPTIDE\d+\{([^}]+)\}/g)].map((match) => match[1]).sort();
}

async function monomerCenters(sk) {
  const centers = [];
  for (let index = 1; index <= ORIGINAL_MONOMERS; index += 1) {
    const rect = await sk.rendered_object_rect('atom', index);
    centers.push({ x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 });
  }
  return centers;
}

test.describe('SKETCH-2911 long monomer labels', () => {
  test('horizontal nucleotide drag leaves a visible gap beside a wide base', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('Choose a custom nucleotide with a wide base label', async () => {
      await clickWidget(page, 'monomeric_btn');
      await clickWidget(page, 'nucleic_monomer_btn');
      await clickWidget(page, 'na_custom_nt_btn');
      await holdPopup(page, 'na_custom_nt_btn');
      await clickWidget(page, 'base_le');
      await page.keyboard.press('ControlOrMeta+a');
      await page.keyboard.type('missingLongBase');
      await expect
        .poll(async () => (await widgetState(page, 'base_le')).text)
        .toBe('missingLongBase');
      await page.keyboard.press('Enter');
      await page.mouse.up();
    });

    const start = await drawingAreaCenter(page);
    await test.step('Drag right to rotate the nucleotide by 90 degrees', async () => {
      // The nucleotide's base starts below its sugar and follows the drag.
      await mouseDrag(page, start, { x: start.x + 100, y: start.y });
    });

    const atoms = await waitAtoms(page, 3);
    expect(atoms.map((atom) => atom.label).sort()).toEqual(['P', 'R', 'missingLongBase']);
    const sugar = atoms.find((atom) => atom.label === 'R');
    const base = atoms.find((atom) => atom.label === 'missingLongBase');

    await test.step('The rendered labels are side by side with a visible gap', async () => {
      expect(
        Math.abs(base.centerY - sugar.centerY),
        'base aligns horizontally with sugar',
      ).toBeLessThanOrEqual(1);
      expect(base.width, 'the long base label is rendered at its full width').toBeGreaterThan(
        sugar.width * 2,
      );
      // Compare measured bead edges so default center spacing cannot hide an
      // overlap. A small pixel margin tolerates integer rounding in the bridge.
      expect(base.x - (sugar.x + sugar.width), 'gap between sugar and base labels').toBeGreaterThan(
        1,
      );
    });

    await test.step('The sugar remains anchored at the mouse press', async () => {
      expect(Math.abs(sugar.centerX - start.x), 'sugar x').toBeLessThanOrEqual(1);
      expect(Math.abs(sugar.centerY - start.y), 'sugar y').toBeLessThanOrEqual(1);
    });
  });

  for (const gesture of ['click', 'drag']) {
    test(`${gesture} adds long-label monomers without moving existing chains`, async ({ page }) => {
      const sk = new Sketcher(page);
      await sk.open();

      await test.step('Import -> Load monomer database -> chembl_monomers_no_id.json', async () => {
        await sk.import_menu('load_monomer_database', MONOMER_DATABASE);
        // Wait for the upload's message dialog, then click its sole OK button.
        // Selecting the ChEMBL-only analog below also verifies the DB loaded.
        const buttons = await widgetRect(page, 'button_box');
        await mouseClick(page, buttons.x + buttons.width - 30, buttons.y + buttons.height / 2);
      });

      await test.step('Draw two separate vertical alanine chains', async () => {
        await clickWidget(page, 'monomeric_btn');
        await clickWidget(page, 'amino_monomer_btn');
        await clickPopupTool(page, 'analog_A_btn');
        const center = await drawingAreaCenter(page);
        for (const [chain, offset] of [-140, 140].entries()) {
          await mouseClick(page, center.x + offset, center.y - 60);
          for (let monomer = 1; monomer < 3; monomer += 1) {
            const rect = await sk.rendered_object_rect('atom', chain * 3 + monomer);
            const start = { x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 };
            // Re-enter the freshly drawn bead so Qt updates its attachment
            // point hover state before the next press.
            await page.mouse.move(start.x + 40, start.y - 40);
            await mouseDrag(page, start, { x: start.x, y: start.y + 60 });
          }
        }
        await expect.poll(() => peptideSequences(page)).toEqual(['A.A.A', 'A.A.A']);
      });

      const before = await monomerCenters(sk);
      // Confirm the setup has two vertical chains before measuring their
      // response to a long label. Derive placement from the rendered chains.
      for (const chain of [before.slice(0, 3), before.slice(3)]) {
        expect(
          Math.max(...chain.map(({ x }) => x)) - Math.min(...chain.map(({ x }) => x)),
        ).toBeLessThan(1);
        expect(Math.abs(chain[2].y - chain[0].y)).toBeGreaterThan(50);
      }
      const middle = { x: (before[1].x + before[4].x) / 2, y: before[1].y };

      await test.step('Long-press a monomer button and choose a long symbol', async () => {
        await clickWidget(page, 'ser_btn');
        await clickPopupTool(page, `analog_${LONG_SYMBOL}_btn`);
        await expect.poll(async () => (await widgetState(page, 'ser_btn')).text).toBe(LONG_SYMBOL);
      });

      await test.step(`${gesture} in the gap adds ${gesture === 'click' ? 'one monomer' : 'two connected monomers'}`, async () => {
        if (gesture === 'click') {
          await mouseClick(page, middle.x, middle.y);
        } else {
          await mouseDrag(page, middle, { x: middle.x, y: middle.y + 60 });
        }
        const addedSequence = Array(gesture === 'click' ? 1 : 2)
          .fill(`[${LONG_SYMBOL}]`)
          .join('.');
        await expect
          .poll(() => peptideSequences(page))
          .toEqual(['A.A.A', 'A.A.A', addedSequence].sort());
      });

      await test.step('Every original monomer remains in its original position', async () => {
        const after = await monomerCenters(sk);
        for (let index = 0; index < ORIGINAL_MONOMERS; index += 1) {
          expect(Math.abs(after[index].x - before[index].x), `monomer ${index + 1} x`).toBeLessThan(
            1,
          );
          expect(Math.abs(after[index].y - before[index].y), `monomer ${index + 1} y`).toBeLessThan(
            1,
          );
        }
      });
    });
  }
});
