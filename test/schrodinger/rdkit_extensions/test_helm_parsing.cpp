
#define BOOST_TEST_MODULE rdkit_extensions_helm_parser

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/helm_parser.h"
#include "schrodinger/rdkit_extensions/helm_examples.h"
#include "test_common.h"

namespace bdata = boost::unit_test::data;

BOOST_DATA_TEST_CASE(TestInvalidHelmInputs, bdata::make(INVALID_EXAMPLES),
                     input_helm)
{
    BOOST_CHECK_THROW(helm::HelmParser(input_helm).parse(),
                      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(TestValidHelmInputs, bdata::make(VALID_EXAMPLES),
                     input_helm)
{
    helm::HelmParser(input_helm).parse();
}

BOOST_DATA_TEST_CASE(TestUnsupportedHelmInputs,
                     bdata::make(UNSUPPORTED_EXAMPLES), input_helm)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(helm::HelmParser(input_helm).parse(),
                                    std::invalid_argument,
                                    "currently unsupported");
}

BOOST_DATA_TEST_CASE(TestBranchMonomerGroups,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{A(C.C)P}$$$$V2.0",
                         "PEPTIDE1{A(C.C)}$$$$V2.0",
                         "PEPTIDE1{A(C.[dC])P}$$$$V2.0",
                     }),
                     input_helm)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(helm::HelmParser(input_helm).parse(),
                                    std::invalid_argument,
                                    "Only one branch monomer is allowed");
}
