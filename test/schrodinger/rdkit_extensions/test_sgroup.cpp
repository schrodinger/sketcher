#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sgroup

#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/sgroup.h"

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string S_GROUP_SMILES = "CCCCCCCCC |Sg:n:3,4,5:1:ht:::|";

// Ensure that get_s_group_atom_bonds returns the expected number of bonds
BOOST_AUTO_TEST_CASE(test_get_s_group_atom_bonds)
{
    auto mol = to_rdkit(S_GROUP_SMILES);
    auto s_group = RDKit::getSubstanceGroups(*mol)[0];
    auto bonds = get_s_group_atom_bonds(s_group);
    BOOST_TEST(bonds.size() == 2);
}

} // namespace rdkit_extensions
} // namespace schrodinger
