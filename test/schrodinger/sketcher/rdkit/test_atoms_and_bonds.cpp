#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE atoms_and_bonds

#include <memory>

#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

#include "../test_common.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::function<RDKit::QueryBond::QUERYBOND_QUERY*()>)

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;

/**
 * Make sure that get_bond_type_and_query_label returns the expected values when
 * passed in queries generated from RDKit QueryOps functions
 */
BOOST_DATA_TEST_CASE(
    test_get_bond_type_and_query_label_from_query_func,
    bdata::make(
        std::vector<std::function<RDKit::QueryBond::QUERYBOND_QUERY*()>>{
            RDKit::makeSingleOrDoubleBondQuery,
            RDKit::makeDoubleOrAromaticBondQuery,
            RDKit::makeSingleOrDoubleOrAromaticBondQuery}) ^
        bdata::make(std::vector<std::string>{"S/D", "D/A", "S/D/A"}),
    query_func, exp_label)
{
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmartsToMol("C=C")};
    auto* bond = mol->getBondWithIdx(0);
    bond->setQuery(query_func());
    auto [bond_type, label] = get_bond_type_and_query_label(bond);
    // since we're going to be labeling this bond, we always want to draw it as
    // a single bond
    BOOST_TEST(bond_type == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(label == exp_label);
}

/**
 * Make sure that get_bond_type_and_query_label returns the expected values when
 * passed in queries generated from SMILES strings
 */
BOOST_DATA_TEST_CASE(
    test_get_bond_type_and_query_label_from_SMARTS,
    bdata::make(std::vector<std::string>{"C=C", "C-,:C", "C@-C", "C:C",
                                         "C!=C"}) ^
        bdata::make(std::vector<std::string>{"", "S/A", "Ring+S", "", "!D"}) ^
        bdata::make(std::vector<RDKit::Bond::BondType>{
            RDKit::Bond::BondType::DOUBLE, RDKit::Bond::BondType::SINGLE,
            RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondType::AROMATIC,
            RDKit::Bond::BondType::SINGLE}),
    smarts, exp_label, exp_bond_type)
{
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmartsToMol(smarts)};
    auto bond = mol->getBondWithIdx(0);
    auto [bond_type, label] = get_bond_type_and_query_label(bond);
    BOOST_TEST(bond_type == exp_bond_type);
    BOOST_TEST(label == exp_label);
}

/**
 * Make sure that get_label_for_bond_query returns "Query" when passed a query
 * that it can't parse
 */
BOOST_AUTO_TEST_CASE(test_get_label_for_bond_query_unrecognized_query)
{
    auto query = RDKit::makeBondInNRingsQuery(3);
    auto label = get_label_for_bond_query(query);
    delete query;
    BOOST_TEST(label == "Query");
}

} // namespace sketcher
} // namespace schrodinger