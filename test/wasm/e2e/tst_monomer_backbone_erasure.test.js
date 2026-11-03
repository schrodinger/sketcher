import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';

// SKETCH-2909 certification steps from the case description.
const PEPTIDE = 'PEPTIDE1{A.A.A.A.A.A}$$$$V2.0';

// The wrapper uses one-based indexes for both bonds and monomers.
const ERASURE_SCENARIOS = [
  {
    target: 'the connection between monomers 3 and 4',
    erase: (sk) => sk.click_bond(3),
    expectedHelm: ['PEPTIDE1{A.A.A}|PEPTIDE2{A.A.A}$$$$V2.0'],
  },
  {
    target: 'monomer 3',
    erase: (sk) => sk.click_atom(3),
    // Chain numbering is immaterial; accept either ordering of the fragments.
    expectedHelm: [
      'PEPTIDE1{A.A}|PEPTIDE2{A.A.A}$$$$V2.0',
      'PEPTIDE1{A.A.A}|PEPTIDE2{A.A}$$$$V2.0',
    ],
  },
];

test.describe('SKETCH-2909 monomer backbone erasure', () => {
  for (const { target, erase, expectedHelm } of ERASURE_SCENARIOS) {
    test(`Copy All As HELM exports two chains after erasing ${target}`, async ({ page }) => {
      const sk = new Sketcher(page);
      await sk.open();

      await test.step('Paste the six-alanine HELM peptide into Sketcher', async () => {
        await sk.import_menu('paste_in_text', PEPTIDE);
      });

      await test.step(`Erase ${target}`, async () => {
        await sk.click_tool('erase');
        await erase(sk);
      });

      await test.step('More Actions -> Copy All As -> HELM copies two separate chains', async () => {
        const helm = await sk.copy_all_as_text('helm');
        expect(expectedHelm).toContain(helm.trim());
      });
    });
  }
});
