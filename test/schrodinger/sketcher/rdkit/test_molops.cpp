#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE molops

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "schrodinger/sketcher/rdkit/molops.h"

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

} // namespace sketcher
} // namespace schrodinger