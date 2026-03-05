
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <optional>
#include <string>

#include <boost/test/unit_test.hpp>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"
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

BOOST_AUTO_TEST_CASE(test_assign_stereochemistry_no_false_chirality_from_3d)
{
    // SKETCH-2706: When a 3D molecule has a carbon with two topologically
    // identical substituents (e.g. the two ring carbons of a cyclopropyl
    // group), assign_stereochemistry must not mark it as a stereocenter
    // even if the 3D geometry produces a non-zero chiral volume.
    //
    // methylcyclopropane (CC1CC1): atom 1 is the junction carbon with
    // neighbors: atom 0 (CH3), atom 2 (ring CH2), atom 3 (ring CH2), and 1
    // implicit H. Atoms 2 and 3 are topologically identical, so atom 1 is
    // NOT a true stereocenter.
    auto mol =
        rdkit_extensions::to_rdkit("CC1CC1", rdkit_extensions::Format::SMILES);

    // Set up a 3D conformer with coordinates that give a non-zero chiral
    // volume at atom 1 (junction carbon). This would cause the legacy stereo
    // perception to incorrectly retain the CHI_TETRAHEDRAL_CW tag.
    auto* conf = new RDKit::Conformer(mol->getNumAtoms());
    conf->set3D(true);
    conf->setAtomPos(0, RDGeom::Point3D(-1.0, -1.0, 0.0)); // CH3
    conf->setAtomPos(1, RDGeom::Point3D(0.0, 0.0, 0.0));   // junction C
    conf->setAtomPos(2, RDGeom::Point3D(1.0, 0.5, 0.5));   // ring CH2 (a)
    conf->setAtomPos(3, RDGeom::Point3D(1.0, -0.5, -0.5)); // ring CH2 (b)
    mol->addConformer(conf, /*assignId=*/true);

    rdkit_extensions::assign_stereochemistry(*mol);

    // The junction carbon must not be marked as a stereocenter: its two
    // cyclopropyl ring neighbors (atoms 2 and 3) have identical CIP ranks.
    auto* junction = mol->getAtomWithIdx(1);
    BOOST_TEST(junction->getChiralTag() ==
               RDKit::Atom::ChiralType::CHI_UNSPECIFIED);
}

} // namespace sketcher
} // namespace schrodinger
