#define BOOST_TEST_MODULE test_biologics_fingerprint

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/rdkit_extensions/biologics_fingerprint.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

using namespace schrodinger::rdkit_extensions::fingerprint;
namespace bdata = boost::unit_test::data;

namespace
{

/// Helper to generate fingerprint from HELM string
auto generate_fp(const std::string& helm_string,
                 const BiologicsFingerprintConfig& config = {})
{
    auto mol = helm::helm_to_rdkit(helm_string);
    return generate_biologics_fingerprint(*mol, config);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(basic_fingerprints)

BOOST_AUTO_TEST_CASE(simple_generation)
{
    auto fp = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0");

    BOOST_CHECK(fp != nullptr);
    BOOST_CHECK_EQUAL(fp->getNumBits(), 2048); // Default size
    BOOST_CHECK_GT(fp->getNumOnBits(), 0);     // Some bits should be set
}

BOOST_AUTO_TEST_CASE(different_sizes)
{
    BiologicsFingerprintConfig config;

    config.fp_size = 512;
    auto fp512 = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0", config);
    BOOST_CHECK_EQUAL(fp512->getNumBits(), 512);

    config.fp_size = 4096;
    auto fp4096 = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0", config);
    BOOST_CHECK_EQUAL(fp4096->getNumBits(), 4096);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(subset_relationships)

// Test data: (query, target, should_be_subset, description)
const std::vector<std::tuple<std::string, std::string, bool, std::string>>
    subset_test_data = {
        {"PEPTIDE1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Single residue subset"},
        {"PEPTIDE1{T}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Missing residue"},
        {"RNA1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different polymer type"},
        {"CHEM1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different polymer type"},
        {"PEPTIDE1{A.C}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Two residue subset (start)"},
        {"PEPTIDE1{C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Two residue subset (middle)"},
        {"PEPTIDE1{A.C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Exact match"},
        {"RNA1{A.C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different biologics type"},
};

BOOST_DATA_TEST_CASE(subset_detection, bdata::make(subset_test_data), query,
                     target, expected, description)
{
    const bool result = is_substructure_fingerprint_match(query, target);

    BOOST_CHECK_MESSAGE(result == expected, description << ": expected "
                                                        << expected << ", got "
                                                        << result);
}

BOOST_AUTO_TEST_CASE(transitive_subsets)
{
    // Test subset chain: A ⊆ A.C ⊆ A.C.G ⊆ A.C.G.T
    const auto helm_a = "PEPTIDE1{A}$$$$V2.0";
    const auto helm_ac = "PEPTIDE1{A.C}$$$$V2.0";
    const auto helm_acg = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto helm_acgt = "PEPTIDE1{A.C.G.T}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(helm_a, helm_ac));
    BOOST_CHECK(is_substructure_fingerprint_match(helm_ac, helm_acg));
    BOOST_CHECK(is_substructure_fingerprint_match(helm_acg, helm_acgt));
    BOOST_CHECK(
        is_substructure_fingerprint_match(helm_a, helm_acgt)); // Transitive
}

BOOST_AUTO_TEST_SUITE_END()
