
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::rdkit_extensions::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::rdkit_extensions::EnhancedStereo>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::nullopt_t)

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Make sure that the stereocenter atom in the provided SMILES has the specified
 * enhanced stereo value
 */
void check_atom_stereochemstry(
    const std::string& smiles,
    const std::optional<EnhancedStereo>& exp_enh_stereo)
{
    auto mol = to_rdkit(smiles, rdkit_extensions::Format::EXTENDED_SMILES);
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
    auto mol = to_rdkit("NC(C)C(=O)O");
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

BOOST_AUTO_TEST_CASE(test_wedgeMolBondse)
{
    // SHARED-11495
    auto mol = to_rdkit("C[C@H](Cl)F", rdkit_extensions::Format::SMILES);
    compute2DCoords(*mol);

    auto b = mol->getBondWithIdx(0);

    // The chirality would be encoded with a dash bond
    wedgeMolBonds(*mol, &mol->getConformer());
    BOOST_REQUIRE_EQUAL(b->getBondDir(), RDKit::Bond::BondDir::BEGINDASH);

    // Check that a misleading mol bond cfg property does
    // not result in the wrong wedging.
    b->setProp(RDKit::common_properties::_MolFileBondCfg, 1); // UP
    wedgeMolBonds(*mol, &mol->getConformer());
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
    auto mol = to_rdkit(cxsmiles, Format::EXTENDED_SMILES);
    BOOST_REQUIRE(mol != nullptr);

    // Generate 2D coordinates
    compute2DCoords(*mol);

    // Set BondDir to UNKNOWN to simulate a wiggly bond
    auto* bond = mol->getBondWithIdx(2);
    BOOST_REQUIRE(bond != nullptr);
    bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);

    // Call wedgeMolBonds - it should preserve the wiggly bond
    wedgeMolBonds(*mol, &mol->getConformer());

    // Verify that BondDir::UNKNOWN was preserved
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);
}

} // namespace rdkit_extensions
} // namespace schrodinger
