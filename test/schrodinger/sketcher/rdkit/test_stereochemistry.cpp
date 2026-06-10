
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::sketcher::EnhancedStereo>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::nullopt_t)

namespace schrodinger
{
namespace sketcher
{

/**
 * Make sure that the stereocenter atom in the provided SMILES has the specified
 * enhanced stereo value
 */
void check_atom_stereochemstry(
    const std::string& smiles,
    const std::optional<EnhancedStereo>& exp_enh_stereo)
{
    auto mol = rdkit_extensions::to_rdkit(
        smiles, rdkit_extensions::Format::EXTENDED_SMILES);
    auto* atom = mol->getAtomWithIdx(1);
    auto enh_stereo = get_enhanced_stereo_for_atom(atom);
    BOOST_TEST(enh_stereo == exp_enh_stereo, smiles);
}

BOOST_AUTO_TEST_CASE(test_get_enhanced_stereo_for_atom)
{
    check_atom_stereochemstry("NC(C)C(=O)O", std::nullopt);
    check_atom_stereochemstry("N[C@H](C)C(=O)O", std::nullopt);
    check_atom_stereochemstry(
        "N[C@H](C)C(=O)O |a:1|",
        EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
    check_atom_stereochemstry(
        "N[C@H](C)C(=O)O |&1:1|",
        EnhancedStereo(RDKit::StereoGroupType::STEREO_AND, 1));
    check_atom_stereochemstry(
        "N[C@H](C)C(=O)O |o2:1|",
        EnhancedStereo(RDKit::StereoGroupType::STEREO_OR, 2));
}

BOOST_AUTO_TEST_CASE(test_set_enhanced_stereo_for_atom)
{
    auto mol = rdkit_extensions::to_rdkit("NC(C)C(=O)O");
    auto* atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) == std::nullopt);
    BOOST_TEST(mol->getStereoGroups().empty());

    // add enhanced stereochemistry for the atom
    set_enhanced_stereo_for_atom(
        atom, EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
    BOOST_TEST(mol->getStereoGroups().size() == 1);

    // switch the stereochemistry for the atom
    set_enhanced_stereo_for_atom(
        atom, EnhancedStereo(RDKit::StereoGroupType::STEREO_AND, 2));
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               EnhancedStereo(RDKit::StereoGroupType::STEREO_AND, 2));
    // make sure that the absolute stereochemistry group was deleted now that
    // it's empty
    BOOST_TEST(mol->getStereoGroups().size() == 1);
}

BOOST_AUTO_TEST_CASE(test_wedgeMolBonds)
{
    // SHARED-11495
    auto mol = rdkit_extensions::to_rdkit("C[C@H](Cl)F",
                                          rdkit_extensions::Format::SMILES);
    rdkit_extensions::compute2DCoords(*mol);

    auto b = mol->getBondWithIdx(0);

    // The chirality would be encoded with a dash bond
    rdkit_extensions::wedgeMolBonds(*mol, &mol->getConformer());
    BOOST_REQUIRE_EQUAL(b->getBondDir(), RDKit::Bond::BondDir::BEGINDASH);

    // Check that a misleading mol bond cfg property does
    // not result in the wrong wedging.
    b->setProp(RDKit::common_properties::_MolFileBondCfg, 1); // UP
    rdkit_extensions::wedgeMolBonds(*mol, &mol->getConformer());
    BOOST_CHECK_EQUAL(b->getBondDir(), RDKit::Bond::BondDir::BEGINDASH);

    // the outdated property should be removed after recalculating wedges
    BOOST_CHECK_EQUAL(b->hasProp(RDKit::common_properties::_MolFileBondCfg),
                      false);
}

/**
 * Verify that wedgeMolBonds preserves wiggly bonds (unknown stereochemistry).
 * Without this preservation, wiggly bonds would be lost when recalculating
 * wedging, since ClearSingleBondDirFlags clears BondDir::UNKNOWN.
 */
BOOST_AUTO_TEST_CASE(test_wedgeMolBonds_preserves_wiggly_bonds)
{
    // Create a molecule with a wiggly bond using CXSMILES
    const std::string cxsmiles = "CCC(C)N |w:2.3|";
    auto mol = rdkit_extensions::to_rdkit(
        cxsmiles, rdkit_extensions::Format::EXTENDED_SMILES);
    BOOST_REQUIRE(mol != nullptr);

    // Generate 2D coordinates
    rdkit_extensions::compute2DCoords(*mol);

    // Set BondDir to UNKNOWN to simulate a wiggly bond
    auto* bond = mol->getBondWithIdx(2);
    BOOST_REQUIRE(bond != nullptr);
    bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);

    // Call wedgeMolBonds - it should preserve the wiggly bond
    rdkit_extensions::wedgeMolBonds(*mol, &mol->getConformer());

    // Verify that BondDir::UNKNOWN was preserved
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);
}

/**
 * Verify that wedgeMolBonds does not steal a wiggly bond to express the stereo
 * of an adjacent chiral atom.
 *
 * WedgeMolBonds prefers terminal-neighbor bonds for wedging. If the wiggly bond
 * connects to a terminal atom (N in this case), WedgeMolBonds would naively
 * steal it to express the adjacent chiral atom's stereo. When we then restore
 * it to UNKNOWN, the chiral atom loses its wedge — the stereo is no longer
 * represented.
 *
 * The fix: before calling WedgeMolBonds, temporarily change wiggly bonds to
 * Bond::OTHER type so WedgeMolBonds skips them entirely (it only considers
 * Bond::SINGLE bonds), then restore them to SINGLE + UNKNOWN afterward.
 */
BOOST_AUTO_TEST_CASE(
    test_wedgeMolBonds_does_not_steal_wiggly_bond_from_adjacent_chiral_atom)
{
    // CC[C@@H](CC)N: atom 2 (C@@H) has defined stereo;
    // the bond to N (atom 5, degree-1 terminal) is the wiggly bond.
    // WedgeMolBonds prefers terminal-neighbor bonds, so without the fix it
    // picks this bond for atom 2's stereo — which we then overwrite with
    // UNKNOWN.
    auto mol = rdkit_extensions::to_rdkit("CC[C@@H](CC)N",
                                          rdkit_extensions::Format::SMILES);
    BOOST_REQUIRE(mol != nullptr);
    rdkit_extensions::compute2DCoords(*mol);

    // Simulate wiggly bond state (as mol_model.cpp does when adding a molecule)
    auto* wiggly_bond = mol->getBondBetweenAtoms(2, 5);
    BOOST_REQUIRE(wiggly_bond != nullptr);
    wiggly_bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);

    rdkit_extensions::wedgeMolBonds(*mol, &mol->getConformer());

    // The wiggly bond must remain UNKNOWN
    BOOST_TEST(wiggly_bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);

    // The chiral atom (C@@H, atom 2) must still have its stereo represented
    // by at least one other bond with a wedge/dash direction.
    auto* stereo_atom = mol->getAtomWithIdx(2);
    bool has_stereo_bond = false;
    for (auto* b : mol->atomBonds(stereo_atom)) {
        if (b != wiggly_bond &&
            (b->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE ||
             b->getBondDir() == RDKit::Bond::BondDir::BEGINDASH)) {
            has_stereo_bond = true;
            break;
        }
    }
    BOOST_TEST(has_stereo_bond,
               "The chiral atom adjacent to the wiggly bond should still have "
               "its stereo expressed by a wedge or dash on another bond");
}

} // namespace sketcher
} // namespace schrodinger
