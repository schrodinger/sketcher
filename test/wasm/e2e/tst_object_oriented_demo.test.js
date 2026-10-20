/**
 * Temporary live demonstration of the Squish-compatible Sketcher wrapper.
 *
 * Aspirin (PubChem CID 2244), canonical SMILES:
 *   CC(=O)Oc1ccccc1C(=O)O
 *
 * This is deliberately an integration demonstration, not a committed
 * regression test.  It composes public wrapper methods only: it does not call
 * a molecule backend API to construct or modify the editor state.
 */
import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';

const ASPIRIN_SMILES = 'CC(=O)Oc1ccccc1C(=O)O';
const DEMO_PAUSE_MS = Number(process.env.PLAYWRIGHT_DEMO_PAUSE_MS || 0);

async function pause(page, label) {
  if (DEMO_PAUSE_MS > 0) {
    // Keep a visible caption in the test trace while allowing a presentation
    // run to linger at each meaningful user-facing state.
    await test.step(`pause: ${label}`, () => page.waitForTimeout(DEMO_PAUSE_MS));
  }
}

function v3000Counts(sdf) {
  const match = sdf.match(/^M  V30 COUNTS (\d+) (\d+) /m);
  return match && { atoms: Number(match[1]), bonds: Number(match[2]) };
}

test('temporary object-oriented wrapper demonstration: aspirin', async ({ page }) => {
  test.setTimeout(120_000);
  const sk = new Sketcher(page);
  await sk.open();

  await test.step('Load aspirin through the visible Import -> Paste in Text workflow', async () => {
    await sk.import_menu('paste_in_text', ASPIRIN_SMILES);
    await pause(page, 'imported aspirin');
  });

  // This uses the existing GUI Copy All As -> SDF path, then parses the
  // resulting editor data into the wrapper's source-compatible structure form.
  const aspirinStructure = await sk.get_structure_information();
  expect(aspirinStructure.atoms).toHaveLength(13);
  expect(aspirinStructure.bonds).toHaveLength(13);

  await test.step('Clear and redraw aspirin with wrapper-managed mouse gestures', async () => {
    await sk.click_button('clear');
    await sk.wait_for_empty_structure();
    await pause(page, 'empty canvas');
    await sk.build_structure(aspirinStructure, { sort: true });
    await pause(page, 'aspirin redrawn from atoms and bonds');
  });

  await test.step('Use live object geometry for additive atom and bond selection', async () => {
    await sk.click_tool('rect_btn');
    for (const index of [1, 2, 3]) await sk.click_atom(index, true, 'shift');
    for (const index of [1, 2]) await sk.click_bond(index, true, 'shift');
    await pause(page, 'geometry-selected atoms and bonds');
  });

  await test.step('Apply and undo a representative More Actions modification', async () => {
    await sk.more_actions_menu('modify_all', 'add_explicit_hydrogens');
    await pause(page, 'explicit hydrogens added');
    await sk.click_button('undo');
  });

  await test.step('Copy an interoperable SDF representation through the visible menu', async () => {
    // The preceding selection step intentionally demonstrates additive object
    // selection. Clear it through the visible Select controls so Copy All As
    // verifies the complete reconstructed molecule rather than that subset.
    await sk.click_button('clear_selection');
    const exportedSdf = await sk.copy_all_as_text('sdf');
    expect(v3000Counts(exportedSdf)).toEqual({ atoms: 13, bonds: 13 });
    expect(exportedSdf).toContain('M  V30 BEGIN ATOM');
    await pause(page, 'final redrawn aspirin');
  });
});
