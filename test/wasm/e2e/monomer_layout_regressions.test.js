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
