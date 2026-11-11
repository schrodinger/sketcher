import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import {
  clear,
  clickItem,
  clickWidget,
  getExportedSmiles,
  getExportedText,
  loadStructure,
  setClipboardText,
  waitAtoms,
} from './e2e_helpers.js';

// SKETCH-2872: wedges and dashes that are already present on an input molecule
// with 2D coordinates must be kept, rather than being cleared and recalculated.
// Copy and paste within the Sketcher round-trips through an RDKit pickle, so it
// relies on this.  The pickle includes any MDL wedging information that the
// molecule still carries, which must not override the wedges that are present.

// Alanine.  The bond indices below are 0-based and follow the SMILES order:
//   0: N-C*   1: C*-CH3   2: C*-C(=O)   3: C=O   4: C-OH
// where C* is the stereocenter.
const ALANINE = 'N[C@H](C)C(=O)O';
const ALANINE_ATOM_COUNT = 6;
const STEREOCENTER_TO_METHYL_BOND = 1;
const STEREOCENTER_TO_CARBOXYL_BOND = 2;
const CARBOXYL_C_TO_OH_BOND = 4;

// Two stereocenters (atoms 2 and 4), each with a wedge that is read from the
// molblock.  The 0-based indices of the bonds around the first stereocenter are:
//   1: C*-F (wedged)   2: C*-C* (the bond to the second stereocenter)
// and the wedge of the second stereocenter is on its C*-Cl bond, which is bond 3.
const TWO_STEREOCENTERS = `
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 F 1.299038 2.250000 0.000000 0
M  V30 4 C 2.598076 0.000000 0.000000 0
M  V30 5 Cl 2.598076 -1.500000 0.000000 0
M  V30 6 C 3.897114 0.750000 0.000000 0
M  V30 7 C 5.196152 0.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 4
M  V30 4 1 4 5 CFG=1
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 2 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
`;
const TWO_STEREOCENTERS_ATOM_COUNT = 7;
const FIRST_STEREOCENTER_TO_F_BOND = 1;
const FIRST_TO_SECOND_STEREOCENTER_BOND = 2;
const SECOND_STEREOCENTER_TO_CL_BOND = 3;

const WEDGE_CFG = 1;
const DASH_CFG = 3;

/**
 * Return every wedged or dashed bond in a V3000 molblock as
 * {bond, begin, end, cfg}, with 0-based bond and atom indices.  The begin atom
 * is included since it determines which end of the bond is narrow.
 */
function stereoBonds(molblock) {
  const lines = molblock.split(/\r?\n/);
  const blockStart = lines.indexOf('M  V30 BEGIN BOND');
  const blockEnd = lines.indexOf('M  V30 END BOND');
  if (blockStart === -1 || blockEnd === -1) {
    // Every structure in this file has bonds, so a missing bond block means
    // that the export failed.  Without this check, we would report that there
    // are no wedges, which could make a comparison pass for the wrong reason.
    throw new Error(`No V3000 bond block found in:\n${molblock}`);
  }
  return lines
    .slice(blockStart + 1, blockEnd)
    .map((line) => line.match(/^M {2}V30 (\d+) \d+ (\d+) (\d+)(?:.* CFG=(\d))?/))
    .filter((match) => match?.[4] !== undefined)
    .map(([, bond, begin, end, cfg]) => ({
      bond: Number(bond) - 1,
      begin: Number(begin) - 1,
      end: Number(end) - 1,
      cfg: Number(cfg),
    }));
}

/** Return the wedge/dash placement of the current Sketcher contents. */
async function currentStereoBonds(page) {
  return stereoBonds(await getExportedText(page, 'MDL_MOLV3000'));
}

/**
 * Return the stereochemistry and the wedge/dash placement of the current
 * Sketcher contents.
 */
async function wedgeState(page) {
  return {
    smiles: await getExportedSmiles(page),
    stereoBonds: await currentStereoBonds(page),
  };
}

/**
 * Wait for the bond with the given index to be drawn with the given CFG value
 * (or to be drawn as a plain bond if cfg is undefined).
 */
async function expectBondCfg(page, bond, cfg) {
  await expect
    .poll(async () => (await currentStereoBonds(page)).find((item) => item.bond === bond)?.cfg, {
      message: `bond ${bond} CFG`,
    })
    .toBe(cfg);
}

/** Add alanine through Import -> Paste in Text */
async function importAlanine(sk, page) {
  await sk.import_menu('paste_in_text', ALANINE);
  await waitAtoms(page, ALANINE_ATOM_COUNT);
}

/**
 * Replace a stereocenter's wedge with a dash on a different bond.
 *
 * @param {import('@playwright/test').Page} page
 * @param {number} wedgedBond - the bond that currently has the wedge
 * @param {number} newBond - the bond to draw the dash on
 */
async function moveStereocenterWedge(page, wedgedBond, newBond) {
  await clickWidget(page, 'stereo_bond2_btn');
  await clickItem(page, 'bond', newBond);
  await expectBondCfg(page, newBond, DASH_CFG);
  await clickWidget(page, 'single_bond_btn');
  await clickItem(page, 'bond', wedgedBond);
  await expectBondCfg(page, wedgedBond, undefined);
}

/**
 * Replace alanine's wedge, which the Sketcher places on the methyl bond, with a
 * dash on the carboxyl bond, which isn't where RDKit would put it.
 */
async function moveAlanineWedge(page) {
  await moveStereocenterWedge(page, STEREOCENTER_TO_METHYL_BOND, STEREOCENTER_TO_CARBOXYL_BOND);
}

/**
 * More Actions -> Copy All, Clear Sketcher, then More Actions -> Paste
 *
 * @param {Sketcher} sk
 * @param {import('@playwright/test').Page} page
 * @param {number} atomCount - the number of atoms in the copied structure
 */
async function copyAllClearAndPaste(sk, page, atomCount = ALANINE_ATOM_COUNT) {
  const sentinel = '__playwright_waiting_for_copy_all__';
  await setClipboardText(page, sentinel);
  await sk.more_actions_menu('copy_all');
  await sk.wait_for_clipboard_change(sentinel, 'Copy All');
  await clear(page);
  await sk.more_actions_menu('paste');
  await waitAtoms(page, atomCount);
}

test.describe('SKETCH-2872: preserving wedges from input molecules', () => {
  let sk;

  test.beforeEach(async ({ page }) => {
    sk = new Sketcher(page);
    await sk.open();
    await importAlanine(sk, page);
  });

  test('a user-added wedge survives Copy All and Paste', async ({ page }) => {
    // These steps match the SKETCH-2872 certification instructions
    await clickWidget(page, 'stereo_bond1_btn');
    await clickItem(page, 'bond', CARBOXYL_C_TO_OH_BOND);
    await expectBondCfg(page, CARBOXYL_C_TO_OH_BOND, WEDGE_CFG);
    const copied = await wedgeState(page);

    await copyAllClearAndPaste(sk, page);
    expect(await wedgeState(page)).toEqual(copied);
  });

  test('a moved stereocenter dash survives Copy All and Paste', async ({ page }) => {
    await moveAlanineWedge(page);
    const copied = await wedgeState(page);
    expect(copied.stereoBonds).toEqual([
      {
        bond: STEREOCENTER_TO_CARBOXYL_BOND,
        begin: 1,
        end: 3,
        cfg: DASH_CFG,
      },
    ]);

    await copyAllClearAndPaste(sk, page);
    expect(await wedgeState(page)).toEqual(copied);
  });

  test('a moved stereocenter dash survives an RDKit pickle import', async ({ page }) => {
    // This is the scenario from the SKETCH-2872 description, where an RDKit
    // molecule with manually placed wedges is passed to the Sketcher
    await moveAlanineWedge(page);
    const exported = await wedgeState(page);
    const pickle = await getExportedText(page, 'RDMOL_BINARY_BASE64');

    await clear(page);
    await loadStructure(page, pickle);
    await waitAtoms(page, ALANINE_ATOM_COUNT);
    expect(await wedgeState(page)).toEqual(exported);
  });

  test('a moved molblock wedge survives Copy All and Paste', async ({ page }) => {
    // A structure that is read from a molblock carries MDL wedging information.
    // Once one of its wedges has been moved, that information is out of date
    // for the edited stereocenter, so it must not be reapplied on paste.
    await clear(page);
    await loadStructure(page, TWO_STEREOCENTERS);
    await waitAtoms(page, TWO_STEREOCENTERS_ATOM_COUNT);
    await expectBondCfg(page, FIRST_STEREOCENTER_TO_F_BOND, WEDGE_CFG);
    await expectBondCfg(page, SECOND_STEREOCENTER_TO_CL_BOND, WEDGE_CFG);

    await moveStereocenterWedge(
      page,
      FIRST_STEREOCENTER_TO_F_BOND,
      FIRST_TO_SECOND_STEREOCENTER_BOND,
    );
    const copied = await wedgeState(page);
    expect(copied.stereoBonds).toEqual([
      { bond: FIRST_TO_SECOND_STEREOCENTER_BOND, begin: 1, end: 3, cfg: DASH_CFG },
      { bond: SECOND_STEREOCENTER_TO_CL_BOND, begin: 3, end: 4, cfg: WEDGE_CFG },
    ]);

    await copyAllClearAndPaste(sk, page, TWO_STEREOCENTERS_ATOM_COUNT);
    expect(await wedgeState(page)).toEqual(copied);
  });
});
