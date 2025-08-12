
#define BOOST_TEST_MODULE test_rdkit_to_helm

#include <boost/algorithm/string.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"
#include "schrodinger/rdkit_extensions/helm_examples.h"

using helm::helm_to_rdkit;
using helm::rdkit_to_helm;

namespace bdata = boost::unit_test::data;

static auto get_roundtripped_helm = [](const auto& input_helm) {
    const auto mol = helm_to_rdkit(input_helm);
    return rdkit_to_helm(*mol);
};

// NOTE: Currently we convert input connections
// with residue names to residue numbers, so those should be skipped
BOOST_DATA_TEST_CASE(TestRoundtrippingValidExamples,
                     bdata::make(VALID_EXAMPLES), input_helm)
{
    // NOTE: This regex is not exhaustive but based on the test examples
    static const std::regex connection_residue_name_regex(R"(\w+\:[R|p])");

    if (std::regex_search(input_helm, connection_residue_name_regex)) {
        return;
    }

    const auto roundtripped_helm = get_roundtripped_helm(input_helm);
    // NOTE: Doing the check this way since some of the inputs do not have the
    // V2.0 suffix
    BOOST_TEST(boost::starts_with(roundtripped_helm, input_helm));
}

// Check that repetitions can be roundtripped
BOOST_DATA_TEST_CASE(TestRoundtrippingMonomerRepetitions,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{A'3'}$$$$V2.0",
                         "PEPTIDE1{A'3-5'}$$$$V2.0",
                         "PEPTIDE1{(X,A)'3'}$$$$V2.0",
                         "PEPTIDE1{(X+A)'3-5'}$$$$V2.0",
                         "PEPTIDE1{(X.A)'3'}$$$$V2.0",
                         "PEPTIDE1{(X.A)'3-5'}$$$$V2.0",
                     }),
                     input_helm)
{
    BOOST_TEST(get_roundtripped_helm(input_helm) == input_helm);
}

// Check that mixtures will be roundtripped
BOOST_DATA_TEST_CASE(TestRoundtrippingMonomerLists,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{(X,A)}$$$$V2.0",
                         "PEPTIDE1{(X+A)}$$$$V2.0",
                         "PEPTIDE1{(X:?,A:0.4)}$$$$V2.0",
                         "PEPTIDE1{(X:?+A:5)}$$$$V2.0",
                     }),
                     input_helm)
{
    BOOST_TEST(get_roundtripped_helm(input_helm) == input_helm);
}

// Check that polymer groups will be roundtripped
BOOST_DATA_TEST_CASE(TestRoundtrippingPolymerGroups,
                     bdata::make(std::vector<std::string>{
                         "CHEM1{*}|CHEM2{*}$$G1(CHEM1,CHEM2)$$V2.0",
                         "CHEM1{*}|CHEM2{*}$$G1(CHEM1:20-30,CHEM2)$$V2.0",
                         "CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0",
                         "CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2:0.3-3)$$V2.0",
                     }),
                     input_helm)
{
    BOOST_TEST(get_roundtripped_helm(input_helm) == input_helm);
}

// Check that connections will be expanded
BOOST_DATA_TEST_CASE(
    TestRoundtrippingConnections,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,5:R2-?:?$$$V2.0",
        "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,A:R3-?:?$$$V2.0",
        "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,1:R3-?:?$$$V2.0",
    }) ^ bdata::make(std::vector<std::string>{
             "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,5:R2-?:?$$$V2.0",
             "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,1:R3-?:?$$$V2.0",
             "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$PEPTIDE1,BLOB1,1:R3-?:?$$$V2.0",
         }),
    input_helm, expected_helm)
{
    BOOST_TEST(get_roundtripped_helm(input_helm) == expected_helm);
}

BOOST_DATA_TEST_CASE(
    TestStrippingSurroundingWhitespace,
    bdata::make(std::vector<std::string>{
        "                 PEPTIDE1{A.A}$$$$V2.0",
        "PEPTIDE1{A.A}$$$$V2.0                ",
        "\n\n\n\n\n\n\n PEPTIDE1{A.A}$$$$V2.0        \n\t      ",
    }),
    input_helm)
{
    std::string expected_helm{"PEPTIDE1{A.A}$$$$V2.0"};
    BOOST_TEST(get_roundtripped_helm(input_helm) == expected_helm);
}
