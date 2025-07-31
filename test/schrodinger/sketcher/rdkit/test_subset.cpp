
#define BOOST_TEST_MODULE subset

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "schrodinger/sketcher/rdkit/subset.h"

using namespace boost::unit_test;

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_get_connected_atoms_and_bonds)
{
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmilesToMol("CCCC.CCC")};
    auto atoms = get_connected_atoms_and_bonds(mol->getAtomWithIdx(0)).first;
    BOOST_CHECK_EQUAL(atoms.size(), 4);
    // Test that the atoms returned by get_connected_atoms_and_bonds are owned
    // by the same molecule as the input atom rather than a copy. SKETCH-2076
    for (auto atom : atoms) {
        BOOST_TEST(&atom->getOwningMol() == mol.get());
    }
}

BOOST_AUTO_TEST_CASE(test_is_contiguous_region)
{
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmilesToMol("CCCC.CC")};
    // create a contiguous region within a disconnected molecule
    std::unordered_set<const RDKit::Atom*> atoms = {
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
        mol->getAtomWithIdx(3)};
    std::unordered_set<const RDKit::Bond*> bonds = {
        mol->getBondBetweenAtoms(0, 1), mol->getBondBetweenAtoms(1, 2),
        mol->getBondBetweenAtoms(2, 3)};
    BOOST_CHECK(is_contiguous_region(atoms, bonds));
    return;

    // remove the first atom from atoms, the region should still be contiguous
    atoms.erase(atoms.begin());
    BOOST_CHECK(is_contiguous_region(atoms, bonds));

    // add it back and remove the first bond from bonds, the region should no
    // longer be contiguous
    atoms.insert(mol->getAtomWithIdx(0));
    bonds.erase(bonds.begin());
    BOOST_CHECK(!is_contiguous_region(atoms, bonds));

    // add all the atoms and bonds in the molecule, the region should not be
    // // contiguous
    atoms = {mol->getAtomWithIdx(0), mol->getAtomWithIdx(1),
             mol->getAtomWithIdx(2), mol->getAtomWithIdx(3),
             mol->getAtomWithIdx(4), mol->getAtomWithIdx(5)};
    bonds = {mol->getBondWithIdx(0), mol->getBondWithIdx(1),
             mol->getBondWithIdx(2), mol->getBondWithIdx(3)};
    BOOST_CHECK(!is_contiguous_region(atoms, bonds));
}

} // namespace sketcher
} // namespace schrodinger