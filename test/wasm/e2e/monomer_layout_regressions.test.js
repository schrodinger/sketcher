import { expect, test } from '@playwright/test';
import {
  clear,
  clickWidget,
  expectLayout,
  getExportedHelm,
  mouseDrag,
  moveAtom,
  pasteText,
  place,
  redo,
  relativeLayout,
  requireRect,
  scene,
  screenshotView,
  undo,
  waitAtoms,
} from './e2e_helpers.js';
import { attachment, connect, openMonomerSketcher } from './monomer_helpers.js';

test.beforeEach(async ({ page }) => {
  await openMonomerSketcher(page);
});
test.setTimeout(120_000);

test('SKETCH-2849: cleanup replaces modestly distorted coordinates and remains stable', async ({
  page,
}) => {
  const initial = 'PEPTIDE1{A.F.C.Y.Q.M.L.C.V.W}$$$$V2.0';
  await pasteText(page, initial);
  await clickWidget(page, 'covalent_or_disulfide_btn');
  await connect(page, 2, 'R3', 7, 'R3');
  await clickWidget(page, 'cleanup_btn');
  const helm = await getExportedHelm(page);
  const expected = await relativeLayout(page);
  // A large distortion made the old implementation discard its coordinates
  // and could mask this bug. Move only a fraction of a bead's dimensions.
  const rect = await requireRect(page, 'atom:9');
  await moveAtom(page, 9, rect.width * 0.6, rect.height * 0.4);
  const distorted = await relativeLayout(page);
  expect(
    Math.hypot(distorted[9].x - expected[9].x, distorted[9].y - expected[9].y),
  ).toBeGreaterThan(1);
  await clickWidget(page, 'cleanup_btn');
  await expectLayout(page, expected);
  expect(await getExportedHelm(page)).toBe(helm);
  await undo(page);
  await expectLayout(page, distorted);
  await redo(page);
  await expectLayout(page, expected);
  await clickWidget(page, 'cleanup_btn');
  await expectLayout(page, expected);
  await clear(page);
  await pasteText(page, helm);
  await expectLayout(page, expected);
});

for (const offset of [0, 1]) {
  test(`SKETCH-2874: joining forward and reverse-drawn chains cleans up like a fresh import (offset ${offset})`, async ({
    page,
  }) => {
    const final = offset
      ? 'PEPTIDE1{G}|PEPTIDE2{C.G.A.A.T.C}$PEPTIDE2,PEPTIDE2,1:R3-6:R3$$$V2.0'
      : 'PEPTIDE1{C.G.A.A.T.C}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$V2.0';
    await pasteText(page, final);
    await clickWidget(page, 'cleanup_btn');
    const expected = await relativeLayout(page);
    await clear(page);
    await pasteText(
      page,
      offset ? 'PEPTIDE1{G}|PEPTIDE2{C.G.A}$$$$V2.0' : 'PEPTIDE1{C.G.A}$$$$V2.0',
    );
    await clickWidget(page, 'cys_btn');
    // Keep the new chain clear of the lower imported polymer when an offset
    // singleton is present. Overlapping beads would target the wrong chain.
    await place(page, -150, 200);
    await waitAtoms(page, offset + 4);
    // Draw C -> N: the insertion order C,T,A differs from sequence A,T,C.
    for (const [button, existing, count] of [
      ['thr_btn', offset + 3, offset + 5],
      ['ala_btn', offset + 4, offset + 6],
    ]) {
      await clickWidget(page, button);
      const start = await attachment(page, existing, 'R1');
      await mouseDrag(page, start, { x: start.x + 110, y: start.y });
      await waitAtoms(page, count);
    }
    await clickWidget(page, 'covalent_or_disulfide_btn');
    await connect(page, offset, 'R3', offset + 3, 'R3');
    await connect(page, offset + 2, 'R2', offset + 5, 'R1');
    expect(await getExportedHelm(page)).toBe(final);
    const order = offset ? [0, 1, 2, 3, 6, 5, 4] : [0, 1, 2, 5, 4, 3];
    const original = await relativeLayout(page, order);
    await clickWidget(page, 'cleanup_btn');
    await expectLayout(page, expected, order);
    expect(await getExportedHelm(page)).toBe(final);
    await undo(page);
    await expectLayout(page, original, order);
    await redo(page);
    await expectLayout(page, expected, order);
    await clickWidget(page, 'cleanup_btn');
    await expectLayout(page, expected, order);
  });
}

test('SKETCH-2822: intrachain connectors avoid bead interiors before and after dragging and image export', async ({
  page,
}) => {
  const sk = new (await import('./sketcher/wrappers/sketcher.js')).Sketcher(page);
  const helm = 'PEPTIDE1{A.K.A.E.G}$PEPTIDE1,PEPTIDE1,2:R3-4:R3$$$V2.0';
  await pasteText(page, helm);
  await clickWidget(page, 'fit_btn');
  await page.mouse.move(1100, 40);
  expect((await scene(page)).connectors).toHaveLength(5);
  for (const connector of (await scene(page)).connectors) {
    expect(
      connector.crossedOtherInteriors,
      `connector ${connector.index} passes through a bead`,
    ).toEqual([]);
  }
  await screenshotView(page, 'intrachain-link.png');
  const before = await relativeLayout(page);
  await moveAtom(page, 3, 35, 65);
  const after = await relativeLayout(page);
  expect(Math.hypot(after[3].x - before[3].x, after[3].y - before[3].y)).toBeGreaterThan(1);
  await page.mouse.move(1100, 40);
  for (const connector of (await scene(page)).connectors) {
    expect(connector.crossedOtherInteriors).toEqual([]);
  }
  expect(await getExportedHelm(page)).toBe(helm);
  await screenshotView(page, 'intrachain-link-after-drag.png');
  const download = await sk.save_image();
  const png = Buffer.from(download.contentBase64, 'base64');
  expect(png.subarray(0, 8)).toEqual(Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]));
  expect(png).toMatchSnapshot('intrachain-link-export.png', { maxDiffPixelRatio: 0.002 });
});

test('SKETCH-2822: an occupied side uses a corner while preserving the backbone', async ({
  page,
}) => {
  // Adjacent residues carry a backbone and a second R3-R3 link. These share a
  // single RDKit bond, but must have two independently rendered connectors.
  await pasteText(page, 'PEPTIDE1{K.E}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0');
  await clickWidget(page, 'fit_btn');
  await page.mouse.move(1100, 40);
  const state = await scene(page);
  expect(state.connectors).toHaveLength(2);
  expect(new Set(state.connectors.map((item) => item.index)).size).toBe(1);
  for (const connector of state.connectors) expect(connector.crossedOtherInteriors).toEqual([]);
  await screenshotView(page, 'adjacent-backbone-and-sidechain.png');
});
