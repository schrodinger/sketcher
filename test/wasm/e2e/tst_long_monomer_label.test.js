import { expect, test } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import {
  choosePopup,
  clickItem,
  clickPopupTool,
  clickWidget as clickE2eWidget,
  expectLayout,
  getExportedHelm,
  holdPopup,
  pasteText,
  place,
  relativeLayout,
  scene,
  screenshotRegion,
  waitAtoms,
} from './e2e_helpers.js';
import { expectDrawingAttachments, loadDatabase, openMonomerSketcher } from './monomer_helpers.js';
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
  test('imported branched peptide has no overlapping monomer labels', async ({ page }) => {
    test.setTimeout(90_000);
    const fixture = JSON.parse(
      await readFile(
        new URL('../../testfiles/long_label_branched_peptide.json', import.meta.url),
        'utf8',
      ),
    );
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('Load the ChEMBL database and import the shared HELM fixture', async () => {
      await loadDatabase(page, 'chembl_monomers_no_id.json');
      await clickWidget(page, 'monomeric_btn');
      await pasteText(page, fixture.helm);
      await waitAtoms(page, fixture.monomer_count);
    });

    const state = await scene(page);
    await test.step('All expected monomers and connections were imported', async () => {
      expect(state.atoms).toHaveLength(fixture.monomer_count);
      expect(state.connectors).toHaveLength(fixture.bond_count);
      const counts = {};
      for (const atom of state.atoms) {
        counts[atom.label] = (counts[atom.label] || 0) + 1;
      }
      expect(counts).toEqual(fixture.monomer_counts);
    });

    await test.step('Every pair of rendered monomer rectangles is separate', async () => {
      const overlaps = [];
      for (let i = 0; i < state.atoms.length; i += 1) {
        const first = state.atoms[i];
        for (const second of state.atoms.slice(i + 1)) {
          const overlapX =
            Math.min(first.x + first.width, second.x + second.width) - Math.max(first.x, second.x);
          const overlapY =
            Math.min(first.y + first.height, second.y + second.height) -
            Math.max(first.y, second.y);
          // The bridge rounds scene rectangles to integer viewport pixels.
          if (overlapX > 1 && overlapY > 1) {
            overlaps.push(
              `${first.index} (${first.label}) and ${second.index} (${second.label}): ` +
                `${overlapX} x ${overlapY} pixel overlap`,
            );
          }
        }
      }
      expect(overlaps, 'overlapping monomer indices and labels').toEqual([]);
    });
  });

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

// Confirm that the predictive hint structure uses the correct connection
// length, even if the existing monomer and predicted monomer both have a long
// name. Make sure that there's enough space between the monomers for the
// connection labels to be visible.
for (const symbol of ['LongAla', 'A']) {
  test(`SKETCH-2871: ${symbol} hover shows attachments and a separated connected preview`, async ({
    page,
  }) => {
    test.setTimeout(90_000);
    await openMonomerSketcher(page);
    if (symbol === 'LongAla') {
      await loadDatabase(page, 'long_symbol.json');
      await choosePopup(page, 'ala_btn', symbol);
    }
    await place(page);
    expect((await waitAtoms(page, 1))[0].label).toBe(symbol);

    // The default C site is labeled on the preview connector, while the free N
    // site remains a nubbin. Inspect the actual painted labels, not the model's
    // attachment definitions. Hovering must not commit the preview molecule.
    await expectDrawingAttachments(page, 0, ['N', 'C'], 'C', ['C', 'N']);
    const hovered = await scene(page);
    expect(hovered.atoms).toHaveLength(1);
    expect(hovered.connectors).toHaveLength(0);
    const bead = hovered.atoms[0];
    const labels = hovered.hintConnections[0].labels;
    for (const label of labels) {
      expect(label.x, `${label.text} is clear of the existing bead`).toBeGreaterThanOrEqual(
        bead.x + bead.width - 1,
      );
    }
    // A tight visual checkpoint also certifies the light hover highlight, the
    // second bead, and the visible connector between the non-overlapping beads.
    await screenshotRegion(page, `${symbol}-hover.png`, {
      x: bead.centerX - 100,
      y: bead.centerY - 60,
      width: 360,
      height: 120,
    });

    await clickItem(page, 'atom', 0);
    const atoms = await waitAtoms(page, 2);
    expect(atoms.map((atom) => atom.label)).toEqual([symbol, symbol]);
    const connected = await scene(page);
    expect(connected.connectors).toHaveLength(1);
    expect(connected.connectors[0]).toMatchObject({ start: 0, end: 1 });
    expect(atoms[1].x - (atoms[0].x + atoms[0].width)).toBeGreaterThan(2);
    const token = symbol === 'A' ? 'A' : `[${symbol}]`;
    expect(await getExportedHelm(page)).toBe(`PEPTIDE1{${token}.${token}}$$$$V2.0`);

    // The committed coordinates must retain their spacing through cleanup.
    const layout = await relativeLayout(page);
    await clickE2eWidget(page, 'cleanup_btn');
    await expectLayout(page, layout);
  });
}
