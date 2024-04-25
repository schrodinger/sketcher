#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_rdkit_to_fasta

#include <boost/test/data/test_case.hpp>
#include "GraphMol/FileParsers/SequenceParsers.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"
#include "schrodinger/rdkit_extensions/fasta/to_string.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"

#include "schrodinger/test/checkexceptionmsg.h" // TEST_CHECK_EXCEPTION_MSG_SUBSTR

namespace bdata = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(TestAtomisticConversion)
{
    auto mol = std::unique_ptr<::RDKit::RWMol>(::RDKit::FASTAToMol("AAPL"));
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(
        fasta::rdkit_to_fasta(*mol), std::invalid_argument,
        "FASTA conversions with atomistic mols are currently unsupported");
}

BOOST_DATA_TEST_CASE(
    TestUnsupportedFeatures,
    bdata::make(std::vector<std::string>{
        R"(PEPTIDE1{A.G.J.K.L}$$$"something"$V2.0)", // extended
                                                     // annotations
        "PEPTIDE1{A.G.J.K.L}|PEPTIDE2{A.G.J.K.L}$$G1(PEPTIDE1+"
        "PEPTIDE2)$$V2.0",                         // polymer groups
        "PEPTIDE1{A.G'3'}$$$$V2.0",                // monomer repeats
        "PEPTIDE1{A.G.J.K.L}|BLOB1{Bead}$$$$V2.0", // blob
                                                   // polymers
        "PEPTIDE1{A.G.J.K.L}|CHEM1{[dR]}$$$$V2.0", // unknown
                                                   // polymers
        "PEPTIDE1{A.G.J.K.L}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2."
        "0", // custom connections
    }),
    input_helm)
{
    auto mol = helm::helm_to_rdkit(input_helm);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(fasta::rdkit_to_fasta(*mol),
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
    BOOST_CHECK_THROW(fasta::rdkit_to_fasta(*mol), std::invalid_argument);
}

BOOST_DATA_TEST_CASE(
    TestValidPeptides,
    bdata::make(std::vector<std::string>{
        ">\nAAZL\n",
        ">something\nARGXCKXEDA\n",
        ">some description\nAAA\n",
        ">some description\nAAA\n>something\nDDD\n",
    }) ^
        bdata::make(std::vector<std::string>{
            "PEPTIDE1{A.A.(E+Q).L}$$$$V2.0",
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
         ">\nAAK\n",
         ">something\nAAB\n",
         ">some description\nAAA\n",
         ">some description\nAAA\n>something\nDDD\n",
     }) ^
     bdata::make(std::vector<std::string>{
         "RNA1{R(A)P.R(A)P.R((G+T+U))P}$$$$V2.0",
         R"(RNA1{R(A)P.R(A)P.R((C+G+T+U))P}"something"$$$$V2.0)",
         R"(RNA1{R(A)P.R(A)P.R(A)P}"some description"$$$$V2.0)",
         R"(RNA1{R(A)P.R(A)P.R(A)P}"some description"|RNA2{R((A+G+T+U))P.R((A+G+T+U))P.R((A+G+T+U))P}"something"$$$$V2.0)",
     })) ^
        bdata::make(std::vector<std::string>{
            "RNA1{[dR](A)P.[dR](A)P.[dR]((G+T+U))P}$$$$V2.0",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR]((C+G+T+U))P}"something"$$$$V2.0)",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR](A)P}"some description"$$$$V2.0)",
            R"(RNA1{[dR](A)P.[dR](A)P.[dR](A)P}"some description"|RNA2{[dR]((A+G+T+U))P.[dR]((A+G+T+U))P.[dR]((A+G+T+U))P}"something"$$$$V2.0)",
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
