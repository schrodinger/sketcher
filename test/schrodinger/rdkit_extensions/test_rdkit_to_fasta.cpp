
#define BOOST_TEST_MODULE test_rdkit_to_fasta

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <rdkit/GraphMol/FileParsers/SequenceParsers.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"
#include "schrodinger/rdkit_extensions/fasta/to_string.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"
#include "schrodinger/test/checkexceptionmsg.h"

namespace bdata = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(TestAtomisticConversion)
{
    auto mol = std::unique_ptr<::RDKit::RWMol>(::RDKit::FASTAToMol("AAPL"));
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(
        std::ignore = fasta::rdkit_to_fasta(*mol), std::invalid_argument,
        "FASTA conversions with atomistic mols are currently unsupported");
}

BOOST_DATA_TEST_CASE(TestUnsupportedFeatures,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{A.G'3'}$$$$V2.0", // monomer repeats
                         "PEPTIDE1{A.G.C.K.L}|BLOB1{Bead}$$$$V2.0", // blob
                                                                    // polymers
                         "PEPTIDE1{A.G.C.K.L}|CHEM1{[dR]}$$$$V2.0", // unknown
                                                                    // polymers
                     }),
                     input_helm)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(std::ignore = fasta::rdkit_to_fasta(*mol),
                                    std::invalid_argument,
                                    "currently unsupported");
}

BOOST_DATA_TEST_CASE(TestUnsupportedNucleotides,
                     bdata::make(std::vector<std::string>{
                         "RNA1{R}$$$$V2.0",     // missing subunit components
                         "RNA1{X(A)P}$$$$V2.0", // unsupported sugar
                         "RNA1{R(A)P.[dR](A)P}$$$$V2.0", // non-uniform sugars
                         "RNA1{R(A)R}$$$$V2.0",          // missing phosphate
                         "RNA1{R.R}$$$$V2.0",            // missing base

                     }),
                     input_helm)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    BOOST_CHECK_THROW(std::ignore = fasta::rdkit_to_fasta(*mol),
                      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(TestUnsupportedMonomers,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{A.A.(E+Q).L}$$$$V2.0",
                         "RNA1{R(A)P.R((C+G+T+U))P}$$$$V2.0",
                         "RNA1{R(Z)P}$$$$V2.0",
                     }),
                     input_helm)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(std::ignore = fasta::rdkit_to_fasta(*mol),
                                    std::invalid_argument,
                                    "Unsupported monomer");
}

BOOST_DATA_TEST_CASE(
    TestValidPeptides,
    bdata::make(std::vector<std::string>{
        ">\nAAOL\n",
        ">something\nARGXCKXEDA\n",
        ">some description\nAAA\n",
        ">some description\nAAA\n>something\nDDD\n",
    }) ^
        bdata::make(std::vector<std::string>{
            "PEPTIDE1{A.A.O.L}$$$$V2.0",
            R"(PEPTIDE1{A.R.G.X.C.K.X.E.D.A}"something"$$$$V2.0)",
            R"(PEPTIDE1{A.A.A}"some description"$$$$V2.0)",
            R"(PEPTIDE1{A.A.A}"some description"|PEPTIDE2{D.D.D}"something"$$$$V2.0)",
        }),
    input_fasta, input_helm)
{
    auto mol = fasta::peptide_fasta_to_rdkit(input_fasta);
    BOOST_TEST(fasta::rdkit_to_fasta(*mol) == input_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == input_helm);
}

BOOST_DATA_TEST_CASE(
    TestValidNucleotides,
    (bdata::make(std::vector<std::string>{
         ">\nAAG\n",
         ">something\nAAU\n",
         ">some description\nAAA\n",
         ">some description\nAAA\n>something\nTTT\n",
     }) ^
     bdata::make(std::vector<std::string>{
         "RNA1{R(A)P.R(A)P.R(G)P}$$$$V2.0",
         R"(RNA1{R(A)P.R(A)P.R(U)P}"something"$$$$V2.0)",
         R"(RNA1{R(A)P.R(A)P.R(A)P}"some description"$$$$V2.0)",
         R"(RNA1{R(A)P.R(A)P.R(A)P}"some description"|RNA2{R(T)P.R(T)P.R(T)P}"something"$$$$V2.0)",
     })) ^
        bdata::make(std::vector<std::string>{
            "RNA1{[dR](A)P.[dR](A)P.[dR](G)P}$$$$V2.0",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR](U)P}"something"$$$$V2.0)",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR](A)P}"some description"$$$$V2.0)",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR](A)P}"some description"|RNA2{[dR](T)P.[dR](T)P.[dR](T)P}"something"$$$$V2.0)",
        }),
    input_fasta, rna_helm, dna_helm)
{
    auto mol = fasta::rna_fasta_to_rdkit(input_fasta);
    BOOST_TEST(fasta::rdkit_to_fasta(*mol) == input_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == rna_helm);

    mol = fasta::dna_fasta_to_rdkit(input_fasta);
    BOOST_TEST(fasta::rdkit_to_fasta(*mol) == input_fasta);
    BOOST_TEST(helm::rdkit_to_helm(*mol) == dna_helm);
}

BOOST_DATA_TEST_CASE(
    TestWritingStructuresWithInterPolymerAndNonlinearConnections,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{A.D(C)P.G}$$$$V2.0",
        "PEPTIDE1{A.G.C.K.L}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2."
        "0",
        R"(PEPTIDE1{A.A.A.C.C}|PEPTIDE2{D.L.L.L.V.V.V}|PEPTIDE3{C.D.D}"Test annotation"$PEPTIDE1,PEPTIDE3,5:R2-1:R1|PEPTIDE3,PEPTIDE2,3:R2-1:R1$$$V2.0)",
    }) ^ bdata::make(std::vector<std::string>{
             ">\nADCPG\n", // branch monomers stripped
             ">\nAGCKL\n", // polymer cycle removed
             ">\nAAACC\n>\nDLLLVVV\n>Test annotation\nCDD\n", // multiple
                                                              // connections
         }),
    input_helm, expected_fasta)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    BOOST_TEST(fasta::rdkit_to_fasta(*mol) == expected_fasta);
}

BOOST_DATA_TEST_CASE(
    TestWritingStructuresWithSupplementaryHELMInformation,
    bdata::make(std::vector<std::string>{
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$G3(PEPTIDE1,PEPTIDE2)$$V2.0)",
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$${"Name":"Test peptides"}$V2.0)",
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$G3(PEPTIDE1,PEPTIDE2)${"Name":"Test peptides"}$V2.0)"}),
    input_helm)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    BOOST_TEST(fasta::rdkit_to_fasta(*mol) == ">\nACD\n>\nDAC\n");
}
