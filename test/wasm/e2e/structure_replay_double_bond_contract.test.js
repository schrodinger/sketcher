import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';

// Acyclovir has several double bonds and exercises the replay state following
// the Single/Double query-bond workflow that originally exposed this defect.
const ACYCLOVIR = 'C1=NC2=C(N1COCCO)NC(=NC2=O)N';
// Aspirin includes the V3000 aromatic bond order (4), which is replayed
// through the visible Aromatic Bond toolbar-popup child.
const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';

function graphSignature(structure) {
  return {
    atoms: structure.atoms.map(({ element, charge }) => ({ element, charge })),
    bonds: structure.bonds.map(({ atom1, atom2, order }) => ({ atom1, atom2, order })),
  };
}

test('structure replay explicitly restores Double Bond after query-bond use', async ({ page }) => {
  test.setTimeout(90_000);
  const sk = new Sketcher(page);
  await sk.open();

  await sk.import_menu('paste_in_text', 'OC1CCCC[C@@H]1Cl');
  await sk.map_imported_atom_indexes();
  await sk.click_tool('single_double', true);
  await sk.click_bond(1);

  await sk.click_button('clear');
  await sk.wait_for_empty_structure();
  await sk.import_menu('paste_in_text', ACYCLOVIR);
  const imported = await sk.get_structure_information();
  await sk.click_button('clear');
  await sk.wait_for_empty_structure();
  sk.replay_tool = null;
  await sk.build_structure(imported);
  expect(graphSignature(await sk.get_structure_information())).toEqual(graphSignature(imported));
});

test('structure replay restores V3000 aromatic bonds through the visible toolbar', async ({
  page,
}) => {
  test.setTimeout(90_000);
  const sk = new Sketcher(page);
  await sk.open();

  await sk.import_menu('paste_in_text', ASPIRIN);
  const imported = await sk.get_structure_information();
  expect(imported.bonds.some(({ order }) => order === 4)).toBe(true);
  await sk.click_button('clear');
  await sk.wait_for_empty_structure();
  sk.replay_tool = null;
  await sk.build_structure(imported);

  expect(graphSignature(await sk.get_structure_information())).toEqual(graphSignature(imported));
});
