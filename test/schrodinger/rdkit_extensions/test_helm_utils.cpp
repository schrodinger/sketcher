
#define BOOST_TEST_MODULE rdkit_extensions_helm_utils

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"

using helm::helm_to_rdkit;
using helm::rdkit_to_helm;
using schrodinger::rdkit_extensions::extract_helm_polymers;
using schrodinger::rdkit_extensions::get_atoms_in_polymer_chain;
using schrodinger::rdkit_extensions::get_atoms_in_polymer_chains;
using schrodinger::rdkit_extensions::get_supplementary_info;
using schrodinger::rdkit_extensions::is_polymer_annotation_s_group;
using schrodinger::rdkit_extensions::is_supplementary_information_s_group;

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
    BOOST_CHECK_THROW(std::ignore = get_atoms_in_polymer_chain(mol, "TEST"),
                      std::invalid_argument);
    BOOST_CHECK_THROW(std::ignore = get_atoms_in_polymer_chains(mol, {}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE(TestPolymerExtraction)

BOOST_AUTO_TEST_CASE(TestAtomisticMol)
{
    ::RDKit::ROMol mol;
    BOOST_CHECK_THROW(std::ignore = extract_helm_polymers(mol, {}),
                      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(TestBadPolymerIds,
                     bdata::make(std::vector<std::string>{"", "UNSUPPORTED"}),
                     polymer_id)
{
    auto mol = helm_to_rdkit(
        "RNA1{[dR](C)P.[dR](A)P}|RNA2{[dR](G)P.[dR](T)P}$$$$V2.0");
    auto extracted_polymers = extract_helm_polymers(*mol, {polymer_id});
    BOOST_TEST(extracted_polymers->getNumAtoms() == 0);
}

BOOST_DATA_TEST_CASE(
    TestExtractionOfPolymers,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{A'3'}$$$$V2.0",
        "PEPTIDE1{G}|PEPTIDE2{A'3'}$$$$V2.0",
        "RNA1{R(C)P}$$$$V2.0",
        "CHEM1{*}$$$$V2.0",
        "BLOB1{Bead}$$$$V2.0",
        "PEPTIDE1{(X:?+A:5)}$$$$V2.0",
        "CHEM1{*}|CHEM2{*}$$$$V2.0",
        "RNA1{[dR](C)P.[dR](A)P}|RNA2{[dR](G)P.[dR](T)P}$$$$V2.0",
        "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$$$$V2.0",
        R"(PEPTIDE1{A'3'}"something here"$$$$V2.0)",
    }) ^ bdata::make(std::vector<std::string>{
             "PEPTIDE1{A'3'}$$$$V2.0",
             "PEPTIDE1{G}$$$$V2.0",
             "RNA1{R(C)P}$$$$V2.0",
             "CHEM1{*}$$$$V2.0",
             "BLOB1{Bead}$$$$V2.0",
             "PEPTIDE1{(X:?+A:5)}$$$$V2.0",
             "CHEM1{*}$$$$V2.0",
             "RNA1{[dR](C)P.[dR](A)P}$$$$V2.0",
             "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$$$$V2.0",
             R"(PEPTIDE1{A'3'}"something here"$$$$V2.0)",
         }),
    input_helm, expected_helm)
{
    auto mol = helm_to_rdkit(input_helm);
    auto extracted_polymers =
        extract_helm_polymers(*mol, {"PEPTIDE1", "RNA1", "CHEM1", "BLOB1"});
    BOOST_TEST(expected_helm == rdkit_to_helm(*extracted_polymers));
}

BOOST_DATA_TEST_CASE(
    TestExtractionOfConnections,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{K.L.C}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0",
        "PEPTIDE1{K.C}|BLOB1{BEAD}$PEPTIDE1,BLOB1,1:R3-?:?$$$V2.0",
        R"(PEPTIDE1{K.C}|BLOB1{BEAD}$PEPTIDE1,BLOB1,1:R3-?:?"Something"$$$V2.0)",
        "PEPTIDE1{K.C}|BLOB2{BEAD}$PEPTIDE1,BLOB2,1:R3-?:?$$$V2.0",
    }) ^
        bdata::make(std::vector<std::string>{
            "PEPTIDE1{K.L.C}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0",
            "PEPTIDE1{K.C}|BLOB1{BEAD}$PEPTIDE1,BLOB1,1:R3-?:?$$$V2.0",
            R"(PEPTIDE1{K.C}|BLOB1{BEAD}$PEPTIDE1,BLOB1,1:R3-?:?"Something"$$$V2.0)",
            "PEPTIDE1{K.C}$$$$V2.0",
        }),
    input_helm, expected_helm)
{
    auto mol = helm_to_rdkit(input_helm);
    auto extracted_polymers =
        extract_helm_polymers(*mol, {"PEPTIDE1", "RNA1", "CHEM1", "BLOB1"});
    BOOST_TEST(expected_helm == rdkit_to_helm(*extracted_polymers));
}

// NOTE: Currently unsupported
BOOST_DATA_TEST_CASE(TestExtractionOfPolymerGroups,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{K.C}|BLOB1{BEAD}$$G1(PEPTIDE1,BLOB1)$$V2.0",
                         "PEPTIDE1{K.C}|BLOB1{BEAD}$$G1(PEPTIDE1+BLOB1)$$V2.0",
                     }),
                     input_helm)
{
    auto mol = helm_to_rdkit(input_helm);
    BOOST_CHECK_THROW(std::ignore = extract_helm_polymers(*mol, {}),
                      std::invalid_argument);
}

// NOTE: Currently unsupported
BOOST_DATA_TEST_CASE(
    TestExtractionOfExtendedAnnotations,
    // clang-format off
    bdata::make(std::vector<std::string>{
        R"(CHEM1{*}|CHEM2{*}$$${"CHEM2":["HELLO"]}$V2.0)",
        R"(CHEM1{*}|CHEM2{*}$$${"CHEM1":["CHEM2"]}$V2.0)",
        R"(CHEM1{*}|CHEM2{*}|CHEM3{*}$$${"KEY":["CHEM2","CHEM3"]}$V2.0)",
        R"(CHEM1{*}|CHEM2{*}$$${")" + ARM_PAIR_KEY + R"(":["CHEM2","CHEM3"]}$V2.0)",
        R"(CHEM1{*}|CHEM2{*}$$${")" + STRAND_PAIR_KEY + R"(":["CHEM2","CHEM3"]}$V2.0)",
        R"(CHEM1{*}|CHEM2{*}$$${"CHEM1": {"CHEM2":"CHEM1"}}$V2.0)",
        }) ^
        bdata::make(std::vector<std::string>{
            R"(CHEM1{*}$$$$V2.0)",
            R"(CHEM1{*}$$$$V2.0)",
            R"(CHEM1{*}$$$$V2.0)",
            R"(CHEM1{*}$$$$V2.0)",
            R"(CHEM1{*}$$$$V2.0)",
            R"(CHEM1{*}$$$$V2.0)",
        }),
    // clang-format on
    input_helm, expected_helm)
{
    auto mol = helm_to_rdkit(input_helm);
    BOOST_CHECK_THROW(std::ignore = extract_helm_polymers(
                          *mol, {"PEPTIDE1", "RNA1", "CHEM1", "BLOB1"}),
                      std::invalid_argument);
}

/**
 * Ensure that we can retrieve and correctly identify a supplementary
 * information S-group
 */
BOOST_AUTO_TEST_CASE(TestSupplementaryInformationSGroup)
{
    std::string helm =
        R"(RNA1{R(A)P.R(C)P.R(G)}$$${"my chain":"my annotation"}$V2.0)";
    auto mol = helm_to_rdkit(helm);
    auto* supplementary_s_group = get_supplementary_info(*mol);
    BOOST_TEST(supplementary_s_group != nullptr);
    BOOST_TEST(is_supplementary_information_s_group(*supplementary_s_group));
    BOOST_TEST(!is_polymer_annotation_s_group(*supplementary_s_group));
}

/**
 * Ensure that we can correctly identify a polymer annotation S-group
 */
BOOST_AUTO_TEST_CASE(TestPolymerAnnotationSGroup)
{
    std::string helm = R"(PEPTIDE1{A.C.D.D.E}"HC"$$$$V2.0)";
    auto mol = helm_to_rdkit(helm);

    // there shouldn't be a supplementary info S-group
    auto* supplementary_s_group = get_supplementary_info(*mol);
    BOOST_TEST(supplementary_s_group == nullptr);

    auto sgroups = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(sgroups.size() == 1);
    auto annotation_s_group = sgroups.front();
    BOOST_TEST(is_polymer_annotation_s_group(annotation_s_group));
    BOOST_TEST(!is_supplementary_information_s_group(annotation_s_group));
}

BOOST_AUTO_TEST_SUITE_END()
