import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';

// SKETCH-2909 certification steps from the case description.
const PEPTIDE = 'PEPTIDE1{A.A.A.A.A.A}$$$$V2.0';
const SPLIT_PEPTIDE = 'PEPTIDE1{A.A.A}|PEPTIDE2{A.A.A}$$$$V2.0';

test.describe('SKETCH-2909 monomer backbone erasure', () => {
  test('Copy All As HELM exports two chains after erasing the middle connection', async ({
    page,
  }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('Paste the six-alanine HELM peptide into Sketcher', async () => {
      await sk.import_menu('paste_in_text', PEPTIDE);
    });

    await test.step('Erase the connection between monomers 3 and 4', async () => {
      await sk.click_tool('erase');
      // The wrapper uses one-based bond indexes; the third backbone bond
      // connects the third and fourth monomers in this linear peptide.
      await sk.click_bond(3);
    });

    await test.step('More Actions -> Copy All As -> HELM copies two separate chains', async () => {
      const helm = await sk.copy_all_as_text('helm');
      expect(helm.trim()).toBe(SPLIT_PEPTIDE);
    });
  });
});
