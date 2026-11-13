
#define BOOST_TEST_MODULE mol_update

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"
#include <algorithm>
#include <optional>
#include <string>

#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/RWMol.h>

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;

std::string LIST_QUERY = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -14.969697 7.242424 0.000000 0
M  V30 2 [C,N] -13.670659 7.992424 0.000000 0
M  V30 3 C -12.371621 7.242424 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$)";

// SKETCH-2710: a molecule with query bonds (V3000 bond type 8 = "any") and a
// stereocenter must not throw when update_molecule_on_change is called. The
// CIPLabeler cannot handle non-integer (UNSPECIFIED) bond orders produced by
// query bonds, so CIP label assignment must be skipped for such molecules.
std::string QUERY_BOND_WITH_STEREO = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 11 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.275976 0.000000 0.000000 0
M  V30 2 C 0.394298 1.213525 0.000000 0
M  V30 3 C -1.032286 0.750000 0.000000 0
M  V30 4 C -1.032286 -0.750000 0.000000 0
M  V30 5 C 0.394298 -1.213525 0.000000 0
M  V30 6 C -2.245812 1.631678 0.000000 0
M  V30 7 C -2.089019 3.123461 0.000000 0
M  V30 8 C -0.718701 3.733566 0.000000 0
M  V30 9 C -3.302545 4.005139 0.000000 0
M  V30 10 N 0.494824 2.851888 0.000000 0
M  V30 11 C -0.561908 5.225349 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 8 1 2
M  V30 2 8 2 3
M  V30 3 8 3 4
M  V30 4 8 4 5
M  V30 5 8 5 1
M  V30 6 1 3 6
M  V30 7 2 6 7
M  V30 8 1 7 8
M  V30 9 1 7 9
M  V30 10 1 8 10
M  V30 11 1 8 11 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 8)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)";

BOOST_AUTO_TEST_CASE(test_update_mol_query_bonds_with_stereo_no_throw)
{
    auto mol = rdkit_extensions::to_rdkit(
        QUERY_BOND_WITH_STEREO, rdkit_extensions::Format::MDL_MOLV3000);
    BOOST_REQUIRE(mol != nullptr);
    // prepare_mol applies wedge bond directions from V3000 CFG fields, which
    // sets CHI tags on atom 3 — needed to trigger CIPLabeler traversal.
    prepare_mol(*mol);
    BOOST_REQUIRE_NO_THROW(update_molecule_on_change(*mol));
}

/**
 * Make sure that prepare_mol flags list queries as dummy atoms
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_list_query)
{
    auto mol = rdkit_extensions::to_rdkit(
        LIST_QUERY, rdkit_extensions::Format::MDL_MOLV3000);
    BOOST_TEST(mol->getNumAtoms() == 3);

    prepare_mol(*mol);

    // the second atom is a list query, check that it's a dummy atom
    BOOST_TEST(mol->getAtomWithIdx(1)->getAtomicNum() ==
               rdkit_extensions::DUMMY_ATOMIC_NUMBER);
}

// A single stereocenter (atom 3) with the wedge on the C3-C4 bond (bond index
// 2)
static const std::string WEDGED_STEREOCENTER = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 2.598076 -0.000000 0.000000 0
M  V30 4 C 3.897114 0.750000 0.000000 0
M  V30 5 O 2.598076 -1.500000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4 CFG=1
M  V30 4 1 3 5
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 3)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)";

// Two stereocenters (atoms 2 and 4, i.e. atom indices 1 and 3). The first has
// its wedge on the C-F bond (bond index 1) and the second has its wedge on the
// C-Cl bond (bond index 3).
static const std::string TWO_WEDGED_STEREOCENTERS = R"(
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
)";

/**
 * @return the number of wedged or dashed bonds in the molecule. If an atom
 * index is given, then only the bonds that start at that atom are counted.
 */
static auto
count_wedges(const RDKit::ROMol& mol,
             const std::optional<unsigned int> begin_atom_idx = std::nullopt)
{
    auto bonds = mol.bonds();
    return std::count_if(bonds.begin(), bonds.end(), [&](auto bond) {
        return rdkit_extensions::is_wedge_or_dash(*bond) &&
               (!begin_atom_idx.has_value() ||
                bond->getBeginAtomIdx() == *begin_atom_idx);
    });
}

/**
 * @return whether the given atom has a tetrahedral chiral tag
 */
static bool is_chiral(const RDKit::Atom* atom)
{
    return atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
           atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW;
}

/**
 * Load WEDGED_STEREOCENTER and remove its wedge, leaving the chiral tag of the
 * stereocenter (atom index 2) in place. No MDL wedging info is left on the
 * molecule.
 */
static auto get_unwedged_mol()
{
    auto mol = rdkit_extensions::to_rdkit(
        WEDGED_STEREOCENTER, rdkit_extensions::Format::MDL_MOLV3000);
    for (auto bond : mol->bonds()) {
        bond->setBondDir(RDKit::Bond::BondDir::NONE);
        bond->clearProp(RDKit::common_properties::_MolFileBondCfg);
    }
    return mol;
}

/**
 * Load WEDGED_STEREOCENTER, then move the wedge from the C3-C4 bond to the C3-O
 * bond (bond index 3) and remove all of the MDL wedging info. (See
 * test_prepare_mol_existing_wedges_override_molblock_wedging for a molecule
 * where the MDL wedging info is left in place.)
 */
static auto get_mol_with_nonstandard_wedge()
{
    auto mol = get_unwedged_mol();
    auto o_bond = mol->getBondWithIdx(3);
    RDKit::Chirality::wedgeBond(o_bond, o_bond->getBeginAtomIdx(),
                                &mol->getConformer());
    return mol;
}

/**
 * Make the given bond a wedge or a dash
 * @param bond The bond to modify
 * @param bond_dir The direction for the bond. Must be BEGINWEDGE or BEGINDASH.
 * @param as_molblock_info If false, the bond direction is set on the bond. If
 * true, the bond direction is cleared and the wedge or dash is instead stored
 * as MDL wedging info, the same as if it had just been read from a V3000
 * molblock.
 */
static void set_wedge(RDKit::Bond* bond, const RDKit::Bond::BondDir bond_dir,
                      const bool as_molblock_info)
{
    if (as_molblock_info) {
        bond->setBondDir(RDKit::Bond::BondDir::NONE);
        int cfg = bond_dir == RDKit::Bond::BondDir::BEGINWEDGE ? 1 : 3;
        bond->setProp(RDKit::common_properties::_MolFileBondCfg, cfg);
    } else {
        bond->setBondDir(bond_dir);
    }
}

/**
 * SKETCH-2872: wedges set on an RDKit mol with 2D coordinates should be
 * honored, even if there's no MDL wedging info present
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_keeps_existing_wedges)
{
    auto mol = get_mol_with_nonstandard_wedge();
    auto o_bond = mol->getBondWithIdx(3);
    auto input_dir = o_bond->getBondDir();
    BOOST_REQUIRE(rdkit_extensions::is_wedge_or_dash(*o_bond));
    auto input_chirality = mol->getAtomWithIdx(2)->getChiralTag();

    prepare_mol(*mol);
    BOOST_TEST(o_bond->getBondDir() == input_dir);
    BOOST_TEST(count_wedges(*mol) == 1);

    // make sure that the stereochemistry is unchanged once the sketcher
    // recalculates it from the wedges
    update_molecule_on_change(*mol);
    BOOST_TEST(o_bond->getBondDir() == input_dir);
    BOOST_TEST(mol->getAtomWithIdx(2)->getChiralTag() == input_chirality);
}

/**
 * SKETCH-2872: If an input mol with 2D coordinates has a wedge that contradicts
 * the chiral tag of the atom that it starts from, then the chiral tag takes
 * precedence. The wedge stays on the same bond, but it's flipped so that it
 * agrees with the chiral tag. Otherwise, a molecule with stale bond directions
 * (e.g. one where the conformer was mirrored after the wedges were assigned)
 * would silently have its stereochemistry inverted. This applies regardless of
 * whether the wedge is stored as a bond direction or as MDL wedging info.
 */
BOOST_DATA_TEST_CASE(test_prepare_mol_chiral_tags_override_conflicting_wedges,
                     bdata::make({false, true}), as_molblock_info)
{
    auto mol = get_mol_with_nonstandard_wedge();
    auto chiral_atom = mol->getAtomWithIdx(2);
    auto input_chirality = chiral_atom->getChiralTag();
    BOOST_REQUIRE(is_chiral(chiral_atom));

    // flip the wedge so that it no longer agrees with the chiral tag
    auto o_bond = mol->getBondWithIdx(3);
    auto correct_dir = o_bond->getBondDir();
    auto conflicting_dir = correct_dir == RDKit::Bond::BondDir::BEGINWEDGE
                               ? RDKit::Bond::BondDir::BEGINDASH
                               : RDKit::Bond::BondDir::BEGINWEDGE;
    set_wedge(o_bond, conflicting_dir, as_molblock_info);

    prepare_mol(*mol);
    BOOST_TEST(o_bond->getBondDir() == correct_dir);
    BOOST_TEST(count_wedges(*mol) == 1);

    update_molecule_on_change(*mol);
    BOOST_TEST(o_bond->getBondDir() == correct_dir);
    BOOST_TEST(chiral_atom->getChiralTag() == input_chirality);
}

/**
 * A wedge that starts at a non-chiral atom and ends at a stereocenter says
 * nothing about that stereocenter. If it's the only wedge touching the
 * stereocenter, then it must not prevent the stereocenter from getting a wedge
 * of its own, since the stereochemistry would be lost otherwise. This applies
 * regardless of whether the wedge is stored as a bond direction or as MDL
 * wedging info.
 */
BOOST_DATA_TEST_CASE(test_prepare_mol_wedge_pointing_into_stereocenter,
                     bdata::make({false, true}), as_molblock_info)
{
    const unsigned int CHIRAL_ATOM_IDX = 2;
    auto mol = get_unwedged_mol();
    auto chiral_atom = mol->getAtomWithIdx(CHIRAL_ATOM_IDX);
    auto input_chirality = chiral_atom->getChiralTag();
    BOOST_REQUIRE(is_chiral(chiral_atom));

    // the C2-C3 bond starts at C2, which isn't a stereocenter
    auto incoming_bond = mol->getBondWithIdx(1);
    BOOST_REQUIRE(incoming_bond->getEndAtomIdx() == CHIRAL_ATOM_IDX);
    set_wedge(incoming_bond, RDKit::Bond::BondDir::BEGINWEDGE,
              as_molblock_info);

    prepare_mol(*mol);
    // the only wedge should now be one that starts at the stereocenter. (Note
    // that RDKit may have picked the C2-C3 bond for that wedge, in which case
    // it will have swapped the bond's start and end atoms.)
    BOOST_TEST(count_wedges(*mol, CHIRAL_ATOM_IDX) == 1);
    BOOST_TEST(count_wedges(*mol) == 1);

    update_molecule_on_change(*mol);
    BOOST_TEST(chiral_atom->getChiralTag() == input_chirality);
}

/**
 * A wedge that starts at a non-chiral atom and ends at a stereocenter should be
 * kept if the stereocenter already has a wedge of its own
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_keeps_extra_wedge_into_stereocenter)
{
    auto mol = get_mol_with_nonstandard_wedge();
    auto chiral_atom = mol->getAtomWithIdx(2);
    auto input_chirality = chiral_atom->getChiralTag();
    auto o_bond = mol->getBondWithIdx(3);
    auto input_dir = o_bond->getBondDir();
    auto incoming_bond = mol->getBondWithIdx(1);
    BOOST_REQUIRE(incoming_bond->getEndAtomIdx() == chiral_atom->getIdx());
    incoming_bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);

    prepare_mol(*mol);
    BOOST_TEST(incoming_bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(incoming_bond->getEndAtomIdx() == chiral_atom->getIdx());
    BOOST_TEST(o_bond->getBondDir() == input_dir);
    BOOST_TEST(count_wedges(*mol) == 2);

    update_molecule_on_change(*mol);
    BOOST_TEST(chiral_atom->getChiralTag() == input_chirality);
}

/**
 * Wavy bonds on an input mol with 2D coordinates should be kept. If the wavy
 * bond was read from a molblock, then its MDL wedging info should be kept as
 * well, since RDKit's CXSMILES writer uses it.
 */
BOOST_DATA_TEST_CASE(test_prepare_mol_keeps_wavy_bonds,
                     bdata::make({false, true}), has_molblock_info)
{
    using RDKit::common_properties::_MolFileBondCfg;
    const int WAVY_CFG = 2;
    auto mol = get_mol_with_nonstandard_wedge();
    auto o_bond = mol->getBondWithIdx(3);
    auto input_dir = o_bond->getBondDir();
    // the C1-C2 bond doesn't touch the stereocenter
    auto wavy_bond = mol->getBondWithIdx(0);
    wavy_bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
    if (has_molblock_info) {
        wavy_bond->setProp(_MolFileBondCfg, WAVY_CFG);
    }

    prepare_mol(*mol);
    BOOST_TEST(wavy_bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);
    BOOST_TEST(o_bond->getBondDir() == input_dir);
    BOOST_TEST(count_wedges(*mol) == 1);
    if (has_molblock_info) {
        int cfg = 0;
        BOOST_TEST(wavy_bond->getPropIfPresent(_MolFileBondCfg, cfg));
        BOOST_TEST(cfg == WAVY_CFG);
    }

    update_molecule_on_change(*mol);
    BOOST_TEST(wavy_bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);
    BOOST_TEST(o_bond->getBondDir() == input_dir);
}

/**
 * If the input has no coordinates, then any existing bond directions are
 * meaningless and should be recalculated
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_recalculates_wedges_without_coords)
{
    const unsigned int CHIRAL_ATOM_IDX = 2;
    auto mol = get_mol_with_nonstandard_wedge();
    auto chiral_atom = mol->getAtomWithIdx(CHIRAL_ATOM_IDX);
    auto input_chirality = chiral_atom->getChiralTag();
    // add a wedge to a bond that doesn't touch a stereocenter, which no
    // wedging calculation would ever produce
    auto non_chiral_bond = mol->getBondWithIdx(0);
    non_chiral_bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);
    mol->clearConformers();

    prepare_mol(*mol);
    BOOST_TEST(non_chiral_bond->getBondDir() == RDKit::Bond::BondDir::NONE);
    // the stereocenter should still have exactly one wedge
    BOOST_TEST(count_wedges(*mol, CHIRAL_ATOM_IDX) == 1);
    BOOST_TEST(count_wedges(*mol) == 1);

    // the new wedge must describe the same stereochemistry using the new
    // coordinates
    update_molecule_on_change(*mol);
    BOOST_TEST(chiral_atom->getChiralTag() == input_chirality);
}

/**
 * SKETCH-2872: A molecule that was read from a molblock still has its MDL
 * wedging info, even if wedges were then set on it using RDKit. Those wedges
 * must take precedence over the MDL wedging info for the stereocenter that they
 * start from, while the MDL wedging info must still be used for the other
 * stereocenters.
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_existing_wedges_override_molblock_wedging)
{
    auto mol = rdkit_extensions::to_rdkit(
        TWO_WEDGED_STEREOCENTERS, rdkit_extensions::Format::MDL_MOLV3000);
    auto first_chiral_atom = mol->getAtomWithIdx(1);
    auto second_chiral_atom = mol->getAtomWithIdx(3);
    BOOST_REQUIRE(is_chiral(first_chiral_atom));
    BOOST_REQUIRE(is_chiral(second_chiral_atom));
    auto first_input_chirality = first_chiral_atom->getChiralTag();
    auto second_input_chirality = second_chiral_atom->getChiralTag();

    auto f_bond = mol->getBondWithIdx(1);
    auto cl_bond = mol->getBondWithIdx(3);
    BOOST_REQUIRE(f_bond->hasProp(RDKit::common_properties::_MolFileBondCfg));
    BOOST_REQUIRE(cl_bond->hasProp(RDKit::common_properties::_MolFileBondCfg));
    // RDKit's molblock parser clears the bond directions once it has
    // determined the chiral tags
    BOOST_REQUIRE(count_wedges(*mol) == 0);

    // put the first stereocenter's wedge on the bond between the two
    // stereocenters, without touching the MDL wedging info
    auto moved_bond = mol->getBondWithIdx(2);
    BOOST_REQUIRE(moved_bond->getBeginAtomIdx() == first_chiral_atom->getIdx());
    RDKit::Chirality::wedgeBond(moved_bond, moved_bond->getBeginAtomIdx(),
                                &mol->getConformer());
    auto moved_dir = moved_bond->getBondDir();
    BOOST_REQUIRE(rdkit_extensions::is_wedge_or_dash(*moved_bond));

    prepare_mol(*mol);
    BOOST_TEST(moved_bond->getBondDir() == moved_dir);
    BOOST_TEST(f_bond->getBondDir() == RDKit::Bond::BondDir::NONE);
    BOOST_TEST(cl_bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(count_wedges(*mol) == 2);

    update_molecule_on_change(*mol);
    BOOST_TEST(moved_bond->getBondDir() == moved_dir);
    BOOST_TEST(count_wedges(*mol) == 2);
    BOOST_TEST(first_chiral_atom->getChiralTag() == first_input_chirality);
    BOOST_TEST(second_chiral_atom->getChiralTag() == second_input_chirality);
}

/**
 * SKETCH-2872: Import a molblock, move one of its wedges, and then round-trip
 * the molecule through an RDKit pickle, which is what happens when copying and
 * pasting within the Sketcher. The pickle includes private properties, so stale
 * MDL wedging info must not be left on the molecule; otherwise it would
 * override the moved wedge when the molecule is pasted.
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_moved_molblock_wedge_survives_pickle)
{
    using RDKit::common_properties::_MolFileBondCfg;
    auto mol = rdkit_extensions::to_rdkit(
        TWO_WEDGED_STEREOCENTERS, rdkit_extensions::Format::MDL_MOLV3000);
    prepare_mol(*mol);
    update_molecule_on_change(*mol);
    auto f_bond = mol->getBondWithIdx(1);
    auto cl_bond = mol->getBondWithIdx(3);
    BOOST_REQUIRE(f_bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_REQUIRE(cl_bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    // the wedges are now stored as bond directions, so the MDL wedging info
    // should be gone
    for (auto bond : mol->bonds()) {
        BOOST_TEST(!bond->hasProp(_MolFileBondCfg));
    }
    auto first_input_chirality = mol->getAtomWithIdx(1)->getChiralTag();
    auto second_input_chirality = mol->getAtomWithIdx(3)->getChiralTag();

    // move the first stereocenter's wedge from the C-F bond to the bond
    // between the two stereocenters. (The Sketcher replaces a bond whenever
    // it's edited, which discards the properties of the original bond.)
    f_bond->setBondDir(RDKit::Bond::BondDir::NONE);
    f_bond->clearProp(_MolFileBondCfg);
    auto moved_bond = mol->getBondWithIdx(2);
    moved_bond->clearProp(_MolFileBondCfg);
    RDKit::Chirality::wedgeBond(moved_bond, moved_bond->getBeginAtomIdx(),
                                &mol->getConformer());
    auto moved_dir = moved_bond->getBondDir();
    BOOST_REQUIRE(rdkit_extensions::is_wedge_or_dash(*moved_bond));
    update_molecule_on_change(*mol);

    auto pickle = rdkit_extensions::to_string(
        *mol, rdkit_extensions::Format::RDMOL_BINARY_BASE64);
    auto pasted_mol = rdkit_extensions::to_rdkit(
        pickle, rdkit_extensions::Format::RDMOL_BINARY_BASE64);
    prepare_mol(*pasted_mol);
    BOOST_TEST(pasted_mol->getBondWithIdx(2)->getBondDir() == moved_dir);
    BOOST_TEST(pasted_mol->getBondWithIdx(1)->getBondDir() ==
               RDKit::Bond::BondDir::NONE);
    BOOST_TEST(pasted_mol->getBondWithIdx(3)->getBondDir() ==
               RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(count_wedges(*pasted_mol) == 2);

    update_molecule_on_change(*pasted_mol);
    BOOST_TEST(count_wedges(*pasted_mol) == 2);
    BOOST_TEST(pasted_mol->getAtomWithIdx(1)->getChiralTag() ==
               first_input_chirality);
    BOOST_TEST(pasted_mol->getAtomWithIdx(3)->getChiralTag() ==
               second_input_chirality);
}

} // namespace sketcher
} // namespace schrodinger
