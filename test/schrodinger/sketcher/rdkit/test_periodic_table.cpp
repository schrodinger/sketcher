
#define BOOST_TEST_MODULE periodic_table

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace boost::unit_test;

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_symbol_to_atomic_number)
{
    std::vector<std::pair<std::string, unsigned int>> symbol_number_pairs = {
        {"H", 1}, {"C", 6},    {"N", 7},  {"O", 8},  {"Og", 118},
        {"o", 8}, {"oG", 118}, {" N", 7}, {"N ", 7},
    };

    for (const auto& symbol_number_pair : symbol_number_pairs) {
        std::string symbol = symbol_number_pair.first;
        unsigned int atomic_number = symbol_number_pair.second;
        BOOST_TEST(symbol_to_atomic_number(symbol) == atomic_number);
        if (atomic_number != 0) {
            BOOST_TEST(
                boost::to_lower_copy(atomic_number_to_symbol(atomic_number)) ==
                boost::trim_copy(boost::to_lower_copy(symbol)));
        }
    }

    // Throws for unknown symbols
    for (const auto& symbol : {"", " ", "unk"}) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(symbol_to_atomic_number(symbol),
                                        Invar::Invariant, "not found");
    }
    // Atom number 0 is a dummy atom in RDKit
    BOOST_TEST(atomic_number_to_symbol(0) == "*");
}

} // namespace sketcher
} // namespace schrodinger