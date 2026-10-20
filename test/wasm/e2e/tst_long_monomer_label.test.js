import { expect, test } from '@playwright/test';
import {
  choosePopup,
  clickItem,
  clickWidget,
  expectLayout,
  getExportedHelm,
  place,
  relativeLayout,
  scene,
  screenshotRegion,
  waitAtoms,
} from './e2e_helpers.js';
import { expectDrawingAttachments, loadDatabase, openMonomerSketcher } from './monomer_helpers.js';

test.setTimeout(90_000);

// Jira's certification says to use a long-name unnatural monomer (e.g. Mono49)
// and compare with a standard amino acid. LongAla is the existing deterministic
// database fixture; select it through the A palette just as Mono49 uses M.
for (const symbol of ['LongAla', 'A']) {
  test(`SKETCH-2871: ${symbol} hover shows attachments and a separated connected preview`, async ({
    page,
  }) => {
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
    await clickWidget(page, 'cleanup_btn');
    await expectLayout(page, layout);
  });
}
