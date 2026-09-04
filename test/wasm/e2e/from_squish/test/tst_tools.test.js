import { expect, test } from '@playwright/test';
import { Sketcher } from '../wrappers/sketcher.js';
import { hideMouseMarker } from '../wrappers/sketcher_wasm.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';
// This is the source-equivalent full visible-tool sweep (roughly 140 real
// popup/canvas interactions), so five minutes is insufficient on Qt/WASM.
test.setTimeout(900_000);

const PERIODIC_TOOLTIPS = [
  ['He', '2 Helium'], ['Li', '3 Lithium'], ['Be', '4 Beryllium'], ['B', '5 Boron'], ['Ne', '10 Neon'],
  ['Na', '11 Sodium'], ['Mg', '12 Magnesium'], ['Al', '13 Aluminium'], ['Si', '14 Silicon'], ['Ar', '18 Argon'],
  ['K', '19 Potassium'], ['Ca', '20 Calcium'], ['Sc', '21 Scandium'], ['Ti', '22 Titanium'], ['V', '23 Vanadium'],
  ['Cr', '24 Chromium'], ['Mn', '25 Manganese'], ['Fe', '26 Iron'], ['Co', '27 Cobalt'], ['Ni', '28 Nickel'],
  ['Cu', '29 Copper'], ['Zn', '30 Zinc'], ['Ga', '31 Gallium'], ['Ge', '32 Germanium'], ['As', '33 Arsenic'],
  ['Se', '34 Selenium'], ['Br', '35 Bromine'], ['Kr', '36 Krypton'], ['Rb', '37 Rubidium'], ['Sr', '38 Strontium'],
  ['Y', '39 Yttrium'], ['Zr', '40 Zirconium'], ['Nb', '41 Niobium'], ['Mo', '42 Molybdenum'], ['Tc', '43 Technetium'],
  ['Ru', '44 Ruthenium'], ['Rh', '45 Rhodium'], ['Pd', '46 Palladium'], ['Ag', '47 Silver'], ['Cd', '48 Cadmium'],
  ['In', '49 Indium'], ['Sn', '50 Tin'], ['Sb', '51 Antimony'], ['Te', '52 Tellurium'], ['I', '53 Iodine'], ['Xe', '54 Xenon'],
  ['Cs', '55 Caesium'], ['Ba', '56 Barium'], ['La', '57 Lanthanum'], ['Ce', '58 Cerium'], ['Pr', '59 Praseodymium'],
  ['Nd', '60 Neodymium'], ['Pm', '61 Promethium'], ['Sm', '62 Samarium'], ['Eu', '63 Europium'], ['Gd', '64 Gadolinium'],
  ['Tb', '65 Terbium'], ['Dy', '66 Dysprosium'], ['Ho', '67 Holmium'], ['Er', '68 Erbium'], ['Tm', '69 Thulium'],
  ['Yb', '70 Ytterbium'], ['Lu', '71 Lutetium'], ['Hf', '72 Hafnium'], ['Ta', '73 Tantalum'], ['W', '74 Tungsten'],
  ['Re', '75 Rhenium'], ['Os', '76 Osmium'], ['Ir', '77 Iridium'], ['Pt', '78 Platinum'], ['Au', '79 Gold'],
  ['Hg', '80 Mercury'], ['Tl', '81 Thallium'], ['Pb', '82 Lead'], ['Bi', '83 Bismuth'], ['Po', '84 Polonium'],
  ['At', '85 Astatine'], ['Rn', '86 Radon'], ['Fr', '87 Francium'], ['Ra', '88 Radium'], ['Ac', '89 Actinium'],
  ['Th', '90 Thorium'], ['Pa', '91 Protactinium'], ['U', '92 Uranium'], ['Np', '93 Neptunium'], ['Pu', '94 Plutonium'],
  ['Am', '95 Americium'], ['Cm', '96 Curium'], ['Bk', '97 Berkelium'], ['Cf', '98 Californium'], ['Es', '99 Einsteinium'],
  ['Fm', '100 Fermium'], ['Md', '101 Mendelevium'], ['No', '102 Nobelium'], ['Lr', '103 Lawrencium'], ['Rf', '104 Rutherfordium'],
  ['Db', '105 Dubnium'], ['Sg', '106 Seaborgium'], ['Bh', '107 Bohrium'], ['Hs', '108 Hassium'], ['Mt', '109 Meitnerium'],
  ['Ds', '110 Darmstadtium'], ['Rg', '111 Roentgenium'], ['Cn', '112 Copernicium'], ['Nh', '113 Nihonium'], ['Fl', '114 Flerovium'],
  ['Mc', '115 Moscovium'], ['Lv', '116 Livermorium'], ['Ts', '117 Tennessine'], ['Og', '118 Oganesson'],
];

async function checkpoint(page, name) {
  await page.mouse.move(0, 0);
  await hideMouseMarker(page);
  await expect(page.locator('#screen canvas')).toHaveScreenshot(`${name}.png`);
}

async function assertToolState(sk, tool, tooltip, clickAndHold = false) {
  const isPeriodicTool = PERIODIC_TOOLTIPS.some(([element]) => element === tool);
  // The source tooltip loop opens/clicks periodic-table children directly,
  // rather than through Sketcher.click_tool(), so it deliberately leaves this
  // wrapper-memory field at its previous value.
  const rememberedPeriodicTool = sk.current_buttons.periodic_table;
  await sk.click_tool(tool, clickAndHold);
  if (isPeriodicTool) sk.current_buttons.periodic_table = rememberedPeriodicTool;
  let state;
  try {
    state = await sk.widget_state(sk.tool_widget_name(tool));
  } catch (error) {
    throw new Error(`Tool state unavailable after selecting ${tool}: ${error.message}`);
  }
  expect(state.enabled).toBe(true);
  expect(state.checked).toBe(true);
  // Periodic-table children close after a browser selection. The persistent
  // last-picked button exposes the same element tooltip without the popup
  // child's atomic-number prefix (for example, "Helium" vs "2 Helium").
  const visibleTooltip = isPeriodicTool
    ? tooltip.replace(/^\d+\s+/, '')
    : tooltip;
  // Bond and selection popup children likewise close after selection; their
  // persistent parent keeps the generic hold-to-change suffix.
  expect(state.toolTip.replace(/ – press & hold to change$/, '')).toBe(
    visibleTooltip.replace(/ – press & hold to change$/, ''),
  );
}

test.describe('tst_tools', () => {
  test('draw-mode basics and atom tools', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source rebuilds only to stabilize desktop Squish coordinates.  The
    // imported V3000 identity map gives the same source atom targets without
    // injecting a separate construction workflow.
    await sk.import_menu('paste_in_text', SOURCE);
    await sk.map_imported_atom_indexes();
    const atom1 = sk.replay_atom_coordinates.get(1);

    // Common element tools increment an existing bond's order.
    await sk.click_tool('C');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_double_C');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_triple_C');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_single_C');

    await sk.click_tool('N');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_double_N');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_triple_N');
    await sk.click_bond(1);
    await checkpoint(page, 'click_bond_to_single_N');

    // A source atom drag adds an attached atom of the selected element.
    await sk.click_tool('C');
    await sk.mouse_drag(atom1.x, atom1.y, -50, 50);
    await checkpoint(page, 'mousedown_drag_C');
    await sk.click_button('undo');
    await sk.click_tool('N');
    await sk.mouse_drag(atom1.x, atom1.y, -50, 50);
    await checkpoint(page, 'mousedown_drag_N');

    await sk.click_button('undo');
    await sk.click_tool('F');
    await sk.click_atom(1);

    // Existing source order: draw triple from atom, then change bond 1.
    await sk.click_tool('triple', true);
    await sk.click_atom(1);
    await checkpoint(page, 'add_bond_of_specified_type');
    await sk.click_tool('triple', true);
    await sk.click_bond(1);
    await checkpoint(page, 'change_bond_to_tool_type');
    await sk.click_tool('coordinate', true);
    await sk.click_bond(1);
    await checkpoint(page, 'change_bond_to_coordinate');
    await sk.click_bond(1);
    await checkpoint(page, 'switch_arrow_coordinate_bond');

    await sk.click_button('clear');
    await sk.click_tool('atom_chain');
    await sk.mouse_drag(-100, 0, 100, 0);
    await checkpoint(page, 'atom_chain');

    // Source's complete visible tool/tooltip loop. It imports a structure
    // first because Select and Move/Rotate are disabled on an empty canvas.
    await sk.import_menu('paste_in_text', SOURCE);
    for (const [tool, tooltip] of [
      ['rect_btn', 'Select structure (with square marquee)'],
      ['lasso_btn', 'Select structure (with lasso)'],
      ['ellipse_btn', 'Select structure (with ellipse marquee)'],
      ['A', 'Any Heavy Atom'], ['AH', 'Any Atom'], ['Q', 'Any Heteroatom'], ['QH', 'Any Heteroatom or H'],
      ['M', 'Any Metal'], ['MH', 'Any Metal or H'], ['X', 'Any Halogen'], ['XH', 'Any Halogen or H'],
    ]) {
      await assertToolState(sk, tool, `${tooltip} – press & hold to change`);
    }
    for (const [tool, tooltip] of PERIODIC_TOOLTIPS) {
      await test.step(`periodic_tool_${tool}`, () => assertToolState(sk, tool, tooltip));
    }
    for (const [tool, tooltip] of [
      ['down', 'Single Down Bond'], ['single_either', 'Single Up or Down Bond'], ['double_either', 'Double Cis or Trans Bond'],
      ['coordinate', 'Coordinate Bond'], ['double', 'Double Bond'], ['triple', 'Triple Bond'], ['zero', 'Zero Order Bond'],
      ['aromatic', 'Aromatic Bond'], ['any', 'Any Bond'], ['single_double', 'Single or Double Bond'],
      ['single_aromatic', 'Single or Aromatic Bond'], ['double_aromatic', 'Double or Aromatic Bond'],
      ['rxn_arrow', 'Reaction Arrow'], ['rxn_plus', 'Reaction Plus'],
      ['add_mapping', 'Map Atoms - Drag from a reactant atom to the matching product atom'],
      ['remove_mapping', 'Delete Mapping - Click mapped atom to remove mapping'],
    ]) {
      await test.step(`popup_tool_${tool}`, () => assertToolState(sk, tool, `${tooltip} – press & hold to change`));
    }
    for (const [tool, tooltip] of [
      ['move_rotate', 'Rotate and translate structure'], ['erase', 'Erase'], ['C', 'Carbon'], ['H', 'Hydrogen'],
      ['N', 'Nitrogen'], ['O', 'Oxygen'], ['P', 'Phosphorus'], ['S', 'Sulphur'], ['F', 'Fluorine'], ['Cl', 'Chlorine'],
      ['cyclohexane', 'Cyclohexane'], ['benzene', 'Benzene'], ['cycloheptane', 'Cycloheptane'], ['cyclopentane', 'Cyclopentane'],
      ['cyclopentadiene', 'Cyclopentadiene'], ['cyclooctane', 'Cyclooctane'], ['cyclobutane', 'Cyclobutane'], ['cyclopropane', 'Cyclopropane'],
      ['up', 'Single Up Bond'], ['single', 'Single Bond'], ['plus_charge', 'Increase Charge'], ['minus_charge', 'Decrease Charge'],
      ['attachment_point', 'Attachment Point'], ['atom_chain', 'Bond Chain - Drag on canvas to draw a chain of carbons'],
    ]) {
      await test.step(`tool_${tool}`, () => assertToolState(sk, tool, tooltip));
    }
    // The source object map registers R-group tooltips last, so this is also
    // the final sticky-tool state Squish asserts below.
    for (const [tool, tooltip] of [
      ['r+', 'New R-Group'], ['r1', 'R-Group 1'], ['r2', 'R-Group 2'], ['r3', 'R-Group 3'], ['r4', 'R-Group 4'],
      ['r5', 'R-Group 5'], ['r6', 'R-Group 6'], ['r7', 'R-Group 7'], ['r8', 'R-Group 8'], ['r9', 'R-Group 9'],
    ]) {
      await test.step(`rgroup_tool_${tool}`, () => assertToolState(sk, tool, `${tooltip} – press & hold to change`));
    }
    expect(sk.current_buttons).toEqual({
      tool: 'r9', wildcard: 'XH', periodic_table: 'Si', stereo: 'double_either',
      bond_order: 'zero', bond_query: 'double_aromatic', 'r-group': 'r9',
      reaction: 'remove_mapping', replace_current_content: 'checked',
      valence_errors: 'checked', heteroatom_colors: 'checked', stereo_labels: 'checked',
      select: 'ellipse_btn',
    });

    await sk.reset_state();
    expect(sk.current_buttons).toEqual({
      tool: 'C', wildcard: 'A', periodic_table: 'Si', stereo: 'down',
      bond_order: 'double', bond_query: 'aromatic', 'r-group': 'r+',
      reaction: 'rxn_arrow', replace_current_content: 'checked',
      valence_errors: 'checked', heteroatom_colors: 'checked', stereo_labels: 'checked',
      select: 'rect_btn',
    });
  });
});
