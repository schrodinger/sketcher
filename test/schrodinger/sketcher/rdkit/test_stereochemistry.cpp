
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

} // namespace sketcher
} // namespace schrodinger
