#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE periodic_table

#include <GraphMol/AtomIterators.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

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

//  All valid atoms
std::vector<std::string> valid_atom_smiles_samples{
    "C",
    "C(C)(C)(C)C",
    "S(C)(C)(C)(C)(C)C",
    "O(C)C",
    "[H+]", // proton
    "[H-]",
    "[HH]",
    "[He]",
    "[He+2]",
    "[C][C] |^5:0,1|",
    "[H][Si] |^5:1|",
    "[CH3+]",
    "[CH3-]",
    "[NH4+]",
    "[Na][H]",
    "[H][Mg][H]",
    "*",             // dummy atom, which also accounts for wildcards
    "*C |$_AP1;$|]", // attachment point
    "[*] |$_R1$|",   // rgroup
    "[Og][Og]([Og])([Og])([Og])([Og])([Og])[Og]",
    "[Lv+4]",
};
// First atom has a valence error!
std::vector<std::string> invalid_first_atom_smiles_samples{
    "[C+5]",
    "C(C)(C)(C)(C)C",
    "S(C)(C)(C)(C)(C)(C)C",
    "[C+](C)(C)(C)C",
    "[C-](C)(C)(C)C",
    "[C](C)(C)(C)C |^1:0|", //  pentavalent due to unpaired electron
    "O(C)=C",
    "[H]",
    "[H+] |^1:0|", // same as [H]
    "[H+] |^2:0|", // non-physical radical count
    "[H+2]",
    "[H-2]",
    "[He+]",
    "[He][He]",
    "[Na]",
    "[Na]([H])[H]",
    "[Mg]",
    "[Mg][H]",
    "[O-3]",
    "[F-2]",
    "[Lv-4]",
};
// Queries never have valence errors
std::vector<std::string> smarts_samples{
    "[#6](C)(C)(C)(C)C",         // query pentavalent carbon
    "[#8](-,=[#6])=[#6]",        // S/D query bond present
    "[!#6&!#1](-[#6])=[#6]",     // Q query atom
    "[#6,#7,#8](-[#6])=[#6]",    // allowed list
    "[!#6&!#7&!#8](-[#6])=[#6]", // disallowed list
    "[#6&R](-[#6])=[#6]",        // advanced query features
};

std::unique_ptr<RDKit::RWMol> to_mol(const std::string& smiles)
{
    int debugParse = 0;
    bool sanitize = false;
    std::unique_ptr<RDKit::RWMol> mol{
        RDKit::SmilesToMol(smiles, debugParse, sanitize)};
    BOOST_REQUIRE(mol != nullptr);
    mol->updatePropertyCache(false);
    return mol;
};

BOOST_DATA_TEST_CASE(test_no_valence_violation,
                     data::make(valid_atom_smiles_samples), smiles)
{
    auto mol = to_mol(smiles);
    for (auto atom : mol->atoms()) {
        BOOST_TEST(!has_valence_violation(atom), smiles);
    }
}

BOOST_DATA_TEST_CASE(test_first_atom_valence_violation,
                     data::make(invalid_first_atom_smiles_samples), smiles)

{
    auto mol = to_mol(smiles);

    // First atom has violation. We don't check the others.
    BOOST_TEST(has_valence_violation(mol->getAtomWithIdx(0)), smiles);
}

BOOST_DATA_TEST_CASE(test_smarts_valence_violation, data::make(smarts_samples),
                     smarts)

{
    std::unique_ptr<RDKit::RWMol> mol{RDKit::SmartsToMol(smarts)};
    BOOST_REQUIRE(mol != nullptr);
    mol->updatePropertyCache(false);
    for (auto atom : mol->atoms()) {
        BOOST_TEST(!has_valence_violation(atom), smarts);
    }
}

} // namespace sketcher
} // namespace schrodinger