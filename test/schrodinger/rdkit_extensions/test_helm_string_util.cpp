#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_test_helm_string_util

#include <boost/test/data/test_case.hpp> // boost::unit_test::data::make
#include <boost/test/unit_test.hpp>

#include <string>

#include "schrodinger/rdkit_extensions/helm/string/util.h"
#include "schrodinger/test/checkexceptionmsg.h" // TEST_CHECK_EXCEPTION_MSG_SUBSTR

BOOST_TEST_DONT_PRINT_LOG_VALUE(HelmVersion)

BOOST_AUTO_TEST_CASE(TestGetHelmVersion)
{
    BOOST_CHECK_EQUAL(get_helm_version("V0.5"), HelmVersion::Unsupported);
    BOOST_CHECK_EQUAL(get_helm_version("$$$$V0.5"), HelmVersion::Unsupported);
    BOOST_CHECK_EQUAL(get_helm_version("$$$$"), HelmVersion::V1);
    BOOST_CHECK_EQUAL(get_helm_version("$$$$V2.0"), HelmVersion::V2);
}

namespace bdata = boost::unit_test::data;

static const std::array<std::string, 4> INVALID_HELM_INPUTS{
    "$$$$V3.0", "$$$$", "POLYMERS$$$$", "POLYMERS}$$$"};
static const std::array<std::string, 4> INVALID_HELM_ERR_MSGS{
    "Unsupported HELM version.", "Input HELM has missing sections.",
    "Input HELM has missing sections.", "Input HELM has missing sections."};

BOOST_DATA_TEST_CASE(TestHelmToHelm2ConversionWithBadInputs,
                     (bdata::make(INVALID_HELM_INPUTS) ^
                      bdata::make(INVALID_HELM_ERR_MSGS)),
                     input_helm, expected_err_msg)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(get_helm2_from_helm(input_helm),
                                    std::invalid_argument, expected_err_msg);
}

static const std::array<std::string, 7> VALID_HELM_INPUTS{
    "POLYMERS}$$$$",
    "POLYMERS}$$$$V2.0",
    "POLYMERS}$CONNECTIONS$$$",
    "POLYMERS}$$H_PAIRINGS$$",
    "POLYMERS}$$$ANNOTATIONS$",
    "POLYMERS}$CONNECTIONS$H_PAIRINGS$$",
    "POLYMERS}$CONNECTIONS$H_PAIRINGS$ANNOTATIONS$"};
static const std::array<std::string, 7> VALID_HELM_OUTPUTS{
    "POLYMERS}$$$$V2.0",
    "POLYMERS}$$$$V2.0",
    "POLYMERS}$CONNECTIONS$$$V2.0",
    "POLYMERS}$H_PAIRINGS$$$V2.0",
    "POLYMERS}$$$ANNOTATIONS$V2.0",
    "POLYMERS}$CONNECTIONS|H_PAIRINGS$$$V2.0",
    "POLYMERS}$CONNECTIONS|H_PAIRINGS$$ANNOTATIONS$V2.0"};

BOOST_DATA_TEST_CASE(TestHelmToHelm2ConversionWithValidInputs,
                     (bdata::make(VALID_HELM_INPUTS) ^
                      bdata::make(VALID_HELM_OUTPUTS)),
                     input_helm, converted_helm)
{
    BOOST_TEST(get_helm2_from_helm(input_helm) == converted_helm);
}
