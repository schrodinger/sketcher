
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
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
 * SKETCH-2759: an MDL/SDF input that explicitly marks a stereogenic double
 * bond as crossed (V3000 CFG=2, V2000 stereo flag 3) must reach the sketcher
 * still tagged STEREOANY/EITHERDOUBLE. update_molecule_on_change re-runs
 * stereo perception with cleanIt=true; without preservation, the explicit
 * "unspecified" intent is lost and the bond is silently assigned E or Z
 * from the 2D layout.
 */
BOOST_AUTO_TEST_CASE(test_paste_preserves_v3000_crossed_double_bond)
{
    static const std::string V3000_CROSSED = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 4.885479 -0.300000 0.000000 0
M  V30 2 C 6.184517 0.450000 0.000000 0
M  V30 3 C 7.483555 -0.300000 0.000000 0
M  V30 4 C 8.782593 0.450000 0.000000 0
M  V30 5 C 10.081631 -0.300000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3 CFG=2
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";

    auto mol = rdkit_extensions::to_rdkit(
        V3000_CROSSED, rdkit_extensions::Format::MDL_MOLV3000);
    update_molecule_on_change(*mol);

    auto* bond = mol->getBondBetweenAtoms(1, 2);
    BOOST_REQUIRE(bond != nullptr);
    BOOST_TEST(bond->getStereo() == RDKit::Bond::BondStereo::STEREOANY);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::EITHERDOUBLE);
}

BOOST_AUTO_TEST_CASE(test_paste_preserves_v2000_crossed_double_bond)
{
    static const std::string V2000_CROSSED = R"CTAB(
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    4.8855   -0.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1845    0.4500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4836   -0.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.7826    0.4500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0816   -0.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  3
  3  4  1  0
  4  5  1  0
M  END
)CTAB";

    auto mol = rdkit_extensions::to_rdkit(
        V2000_CROSSED, rdkit_extensions::Format::MDL_MOLV2000);
    update_molecule_on_change(*mol);

    auto* bond = mol->getBondBetweenAtoms(1, 2);
    BOOST_REQUIRE(bond != nullptr);
    BOOST_TEST(bond->getStereo() == RDKit::Bond::BondStereo::STEREOANY);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::EITHERDOUBLE);
}

BOOST_AUTO_TEST_CASE(test_smiles_specified_stereo_not_clobbered)
{
    auto mol =
        rdkit_extensions::to_rdkit("C/C=C/C", rdkit_extensions::Format::SMILES);
    rdkit_extensions::compute2DCoords(*mol);
    update_molecule_on_change(*mol);

    auto* bond = mol->getBondBetweenAtoms(1, 2);
    BOOST_REQUIRE(bond != nullptr);
    BOOST_TEST(bond->getStereo() != RDKit::Bond::BondStereo::STEREOANY);
    BOOST_TEST(bond->getBondDir() != RDKit::Bond::BondDir::EITHERDOUBLE);
}

} // namespace sketcher
} // namespace schrodinger
