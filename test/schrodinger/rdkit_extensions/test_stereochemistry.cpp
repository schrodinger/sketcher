#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
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

} // namespace rdkit_extensions
} // namespace schrodinger
