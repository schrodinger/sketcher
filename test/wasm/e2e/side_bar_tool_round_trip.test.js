import { expect, test } from '@playwright/test';
import { clickWidget, place, waitAtoms, widgetState } from './e2e_helpers.js';
import { BASE, CHEM, PEPTIDE, defineCustom, openMonomerSketcher } from './monomer_helpers.js';

// Switching the side bar to the other monomer page, or to the atomistic tools,
// and straight back should leave the tool selection exactly as it was. The
// connection and custom monomer tools used to be replaced by the page's
// standard monomer tool, and so did the shared select and erase tools on a
// page switch.

test.beforeEach(async ({ page }) => {
  await openMonomerSketcher(page);
});
test.setTimeout(90_000);

async function expectChecked(page, name, checked) {
  await expect.poll(async () => (await widgetState(page, name)).checked).toBe(checked);
}

test('a page round trip restores the connection tool', async ({ page }) => {
  await clickWidget(page, 'hbond_btn');
  await expectChecked(page, 'hbond_btn', true);
  await clickWidget(page, 'nucleic_monomer_btn');
  await expectChecked(page, 'hbond_btn', false);
  await clickWidget(page, 'amino_monomer_btn');
  await expectChecked(page, 'hbond_btn', true);
});

for (const sharedTool of ['select_tool_btn', 'erase_btn']) {
  test(`a page round trip keeps the shared ${sharedTool} tool active`, async ({ page }) => {
    // The select and erase tools are only enabled when the scene has contents.
    // Ala is already the active tool. (Clicking its checked button would open
    // its analog popup instead, and the canvas click would only close it.)
    await expectChecked(page, 'ala_btn', true);
    await place(page);
    await waitAtoms(page, 1);
    await clickWidget(page, sharedTool);
    await expectChecked(page, sharedTool, true);
    await clickWidget(page, 'nucleic_monomer_btn');
    await expectChecked(page, sharedTool, true);
    await clickWidget(page, 'amino_monomer_btn');
    await expectChecked(page, sharedTool, true);
    await expectChecked(page, 'ala_btn', false);

    // the page's own monomeric tool is still remembered
    await clickWidget(page, 'nucleic_monomer_btn');
    await clickWidget(page, 'na_c_btn');
    await clickWidget(page, 'amino_monomer_btn');
    await expectChecked(page, 'ala_btn', true);
    await expectChecked(page, sharedTool, false);
  });
}

test('a page round trip restores a custom monomer unless the other page replaced it', async ({
  page,
}) => {
  await defineCustom(page, 'peptide', PEPTIDE);
  await expectChecked(page, 'custom_monomer_btn', true);
  await clickWidget(page, 'nucleic_monomer_btn');
  await expectChecked(page, 'custom_monomer_btn', false);
  await clickWidget(page, 'amino_monomer_btn');
  await expectChecked(page, 'custom_monomer_btn', true);
  await place(page);
  await waitAtoms(page, 1);

  // A custom base replaces the stored custom peptide, so the amino acid page
  // falls back to its standard tool, while the nucleic acid page keeps the base
  await defineCustom(page, 'base', BASE);
  await expectChecked(page, 'custom_monomer_btn', true);
  await clickWidget(page, 'amino_monomer_btn');
  await expectChecked(page, 'custom_monomer_btn', false);
  await expectChecked(page, 'ala_btn', true);
  await clickWidget(page, 'nucleic_monomer_btn');
  await expectChecked(page, 'custom_monomer_btn', true);
});

test('a round trip through the atomistic tools restores the connection and custom monomer tools', async ({
  page,
}) => {
  await clickWidget(page, 'covalent_or_disulfide_btn');
  await expectChecked(page, 'covalent_or_disulfide_btn', true);
  await clickWidget(page, 'atomistic_btn');
  await clickWidget(page, 'monomeric_btn');
  await expectChecked(page, 'covalent_or_disulfide_btn', true);

  await defineCustom(page, 'chem', CHEM);
  await expectChecked(page, 'custom_monomer_btn', true);
  await clickWidget(page, 'atomistic_btn');
  await clickWidget(page, 'monomeric_btn');
  await expectChecked(page, 'custom_monomer_btn', true);
  await expectChecked(page, 'covalent_or_disulfide_btn', false);
  await place(page);
  await waitAtoms(page, 1);
});
