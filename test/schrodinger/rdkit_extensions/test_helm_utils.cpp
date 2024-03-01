#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_helm_utils

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

using schrodinger::rdkit_extensions::get_atoms_in_polymer_chain;
using schrodinger::rdkit_extensions::get_atoms_in_polymer_chains;

namespace bdata = boost::unit_test::data;

using atoms_t = std::vector<unsigned int>;
BOOST_TEST_DONT_PRINT_LOG_VALUE(atoms_t)

BOOST_DATA_TEST_CASE(TestGetAtomsInPolymerChainWithSingleId,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1",
                         "PEPTIDE2",
                         "PEPTIDE3",
                     }) ^ bdata::make(std::vector<atoms_t>{
                              {0, 1, 2},
                              {3, 4, 5},
                              {},
                          }),
                     polymer_id, expected_atoms)
{
    std::string input_helm{"PEPTIDE1{A.A.A}|PEPTIDE2{C.C.C}$$$$V2.0"};
    auto mol = helm::helm_to_rdkit(input_helm);
    auto extracted_atoms = get_atoms_in_polymer_chain(*mol, polymer_id);
    BOOST_TEST(extracted_atoms == expected_atoms);
}

BOOST_TEST_DONT_PRINT_LOG_VALUE(std::vector<std::string_view>)
BOOST_DATA_TEST_CASE(TestGetAtomsInPolymerChainWithMultipleIds,
                     bdata::make(std::vector<std::vector<std::string_view>>{
                         {"PEPTIDE1"},
                         {"PEPTIDE1", "PEPTIDE1"},
                         {"PEPTIDE2"},
                         {"PEPTIDE1", "PEPTIDE2"},
                         {"PEPTIDE1", "PEPTIDE3"},
                         {"PEPTIDE3"},
                         {},
                     }) ^ bdata::make(std::vector<atoms_t>{
                              {0, 1, 2},
                              {0, 1, 2},
                              {3, 4, 5},
                              {0, 1, 2, 3, 4, 5},
                              {0, 1, 2},
                              {},
                              {},
                          }),
                     polymer_ids, expected_atoms)
{
    std::string input_helm{"PEPTIDE1{A.A.A}|PEPTIDE2{C.C.C}$$$$V2.0"};
    auto mol = helm::helm_to_rdkit(input_helm);
    auto extracted_atoms = get_atoms_in_polymer_chains(*mol, polymer_ids);
    BOOST_TEST(extracted_atoms == expected_atoms);
}

BOOST_AUTO_TEST_CASE(TestAtomisticMolsAreUnsupported)
{
    ::RDKit::ROMol mol;
    BOOST_CHECK_THROW(get_atoms_in_polymer_chain(mol, "TEST"),
                      std::invalid_argument);
    BOOST_CHECK_THROW(get_atoms_in_polymer_chains(mol, {}),
                      std::invalid_argument);
}
