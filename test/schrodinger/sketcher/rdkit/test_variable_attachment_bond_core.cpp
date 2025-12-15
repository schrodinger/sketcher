
#define BOOST_TEST_MODULE test_variable_attachment_bond

#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"

namespace schrodinger
{
namespace sketcher
{

const std::string VAR_ATTACH_MOL = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.065751 0.028571 0.000000 0
M  V30 2 C 1.302930 0.742857 0.000000 0
M  V30 3 C 1.302930 2.171429 0.000000 0
M  V30 4 C 0.065751 2.885714 0.000000 0
M  V30 5 C -1.171429 2.171429 0.000000 0
M  V30 6 C -1.171429 0.742857 0.000000 0
M  V30 7 C 1.398626 3.765751 0.000000 0
M  V30 8 * 0.327197 1.909982 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 1 8 7 ENDPTS=(3 2 3 4) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$

)";

/**
 * Ensure that add_variable_attachment_bond_to_mol() adds a correctly formed
 * bond to the molecule.  Also make sure that the resulting bond returns the
 * expected values from is_variable_attachment_bond and
 * get_variable_attachment_atoms.
 */
BOOST_AUTO_TEST_CASE(test_add_variable_attachment_bond_to_mol)
{
    auto mol = rdkit_extensions::to_rdkit("C1CCCCC1");
    rdkit_extensions::compute2DCoords(*mol);
    const std::unordered_set<const RDKit::Atom*> atoms{
        mol->getAtomWithIdx(1), mol->getAtomWithIdx(2), mol->getAtomWithIdx(3)};
    auto [dummy_atom, carbon_atom, bond] =
        add_variable_attachment_bond_to_mol(*mol, atoms);
    BOOST_TEST(dummy_atom->getAtomicNum() == 0);
    BOOST_TEST(carbon_atom->getAtomicNum() == 6);
    BOOST_TEST(bond->hasProp(RDKit::common_properties::_MolFileBondAttach));
    BOOST_TEST(bond->hasProp(RDKit::common_properties::_MolFileBondEndPts));
    BOOST_TEST(bond->getProp<std::string>(
                   RDKit::common_properties::_MolFileBondAttach) == "ANY");
    // molfile indices are 1-based instead of 0-based, which is why these
    // numbers are different than the indices above
    BOOST_TEST(bond->getProp<std::string>(
                   RDKit::common_properties::_MolFileBondEndPts) ==
               "(3 2 3 4)");
    BOOST_TEST(bond->getOtherAtom(dummy_atom) == carbon_atom);
    BOOST_TEST(bond->getOtherAtom(carbon_atom) == dummy_atom);

    BOOST_TEST(is_variable_attachment_bond(bond));
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);
    BOOST_TEST(is_dummy_atom_for_variable_attachment_bond(dummy_atom));
    BOOST_TEST(!is_dummy_atom_for_variable_attachment_bond(carbon_atom));
    BOOST_TEST(
        !is_dummy_atom_for_variable_attachment_bond(mol->getAtomWithIdx(1)));
}

/**
 * Ensure that add_variable_attachment_bond_to_mol() adds a correctly formed
 * bond to the molecule, like the test above, but include the atoms with the
 * first and last index to make sure that we don't have any off-by-one errors.
 */
BOOST_AUTO_TEST_CASE(
    test_add_variable_attachment_bond_to_mol_using_first_and_last_atoms)
{
    auto mol = rdkit_extensions::to_rdkit("C1CCCCC1");
    rdkit_extensions::compute2DCoords(*mol);
    const std::unordered_set<const RDKit::Atom*> atoms{
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(3), mol->getAtomWithIdx(5)};
    auto [dummy_atom, carbon_atom, bond] =
        add_variable_attachment_bond_to_mol(*mol, atoms);
    BOOST_TEST(dummy_atom->getAtomicNum() == 0);
    BOOST_TEST(carbon_atom->getAtomicNum() == 6);
    BOOST_TEST(bond->hasProp(RDKit::common_properties::_MolFileBondAttach));
    BOOST_TEST(bond->hasProp(RDKit::common_properties::_MolFileBondEndPts));
    BOOST_TEST(bond->getProp<std::string>(
                   RDKit::common_properties::_MolFileBondAttach) == "ANY");
    // molfile indices are 1-based instead of 0-based, which is why these
    // numbers are different than the indices above
    BOOST_TEST(bond->getProp<std::string>(
                   RDKit::common_properties::_MolFileBondEndPts) ==
               "(3 1 4 6)");
    BOOST_TEST(bond->getOtherAtom(dummy_atom) == carbon_atom);
    BOOST_TEST(bond->getOtherAtom(carbon_atom) == dummy_atom);

    BOOST_TEST(is_variable_attachment_bond(bond));
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);
}

/**
 * Ensure that is_variable_attachment_bond and get_variable_attachment_atoms
 * return the expected results using bonds from an RDKit-created molecule.
 */
BOOST_AUTO_TEST_CASE(test_get_variable_attachment_atoms)
{
    auto mol = rdkit_extensions::to_rdkit(VAR_ATTACH_MOL);
    const auto* normal_bond = mol->getBondWithIdx(0);
    BOOST_TEST(!is_variable_attachment_bond(normal_bond));
    BOOST_TEST(get_variable_attachment_atoms(normal_bond).empty());
    const auto* variable_attachment_bond = mol->getBondWithIdx(6);
    BOOST_TEST(is_variable_attachment_bond(variable_attachment_bond));
    auto atoms = get_variable_attachment_atoms(variable_attachment_bond);
    BOOST_TEST(atoms.size() == 3);
    const std::unordered_set<const RDKit::Atom*> exp_atoms{
        mol->getAtomWithIdx(1), mol->getAtomWithIdx(2), mol->getAtomWithIdx(3)};
    BOOST_TEST(atoms == exp_atoms);
}

/**
 * Ensure that get_variable_attachment_atoms can parse ENDPTS property values
 * with incorrect spaces, and make sure that it doesn't crash if the ENDPTS
 * property can't be parsed
 */
BOOST_AUTO_TEST_CASE(test_get_variable_attachment_atoms_misformatted_bond_props)
{
    auto mol = rdkit_extensions::to_rdkit("C1CCCCC1");
    rdkit_extensions::compute2DCoords(*mol);
    const std::unordered_set<const RDKit::Atom*> atoms{
        mol->getAtomWithIdx(1), mol->getAtomWithIdx(2), mol->getAtomWithIdx(3)};
    auto [dummy_atom, carbon_atom, bond] =
        add_variable_attachment_bond_to_mol(*mol, atoms);
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);

    // extra whitespace should be parseable
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts,
                  "( 3 2    3 4  )");
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);

    // any other errors shouldn't be parseable, but
    // get_variable_attachment_atoms shouldn't crash

    // negative number of atoms
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(-3 2 3 4)");
    const std::unordered_set<const RDKit::Atom*> empty_set;
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    // negative atom index
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(3 2 -3 4)");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    // letter
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(3 2a 3 4)");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    // invalid atom indices
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(3 2 3 17)");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(3 0 3 4 )");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "(3 2 3 9)");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    // missing parenthesis
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts, "3 2 3 4");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
    // text
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts,
                  "This is wrong");
    BOOST_TEST(get_variable_attachment_atoms(bond) == empty_set);
}

/**
 * Ensure that get_variable_attachment_atoms can correctly parse ENDPTS property
 * values containing double-digit numbers.
 */
BOOST_AUTO_TEST_CASE(test_get_variable_attachment_atoms_double_digit_numbers)
{
    auto mol = rdkit_extensions::to_rdkit("CCCCCCCCCCC1CCCCC1");
    rdkit_extensions::compute2DCoords(*mol);
    const std::unordered_set<const RDKit::Atom*> atoms{mol->getAtomWithIdx(11),
                                                       mol->getAtomWithIdx(12),
                                                       mol->getAtomWithIdx(13)};
    auto [dummy_atom, carbon_atom, bond] =
        add_variable_attachment_bond_to_mol(*mol, atoms);
    BOOST_TEST(bond->getProp<std::string>(
                   RDKit::common_properties::_MolFileBondEndPts) ==
               "(3 12 13 14)");
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);
}

} // namespace sketcher
} // namespace schrodinger
