
#define BOOST_TEST_MODULE mol_update

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/constants.h"
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

} // namespace sketcher
} // namespace schrodinger