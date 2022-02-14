/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: text block <-> rdkit mol
 conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_convert

#include <map>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

const std::array<Format, 7> TEXT_FORMATS = {
    Format::SMILES,       Format::EXTENDED_SMILES, Format::SMARTS,
    Format::MDL_MOLV2000, Format::MDL_MOLV3000,    Format::INCHI,
    Format::PDB,
};

const std::array<Format, 4> REACTION_TEXT_FORMATS = {
    Format::SMILES, Format::SMARTS, Format::MDL_MOLV2000, Format::MDL_MOLV3000};

BOOST_AUTO_TEST_CASE(temp)
{
    BOOST_TEST(true);
}

BOOST_DATA_TEST_CASE(test_auto_detect,
                     boost::unit_test::data::make(TEXT_FORMATS))
{
    auto mol = std::shared_ptr<RDKit::RWMol>(RDKit::SmilesToMol("c1ccccc1"));
    auto text = rdmol_to_text(*mol, sample);

    // Check roundtripping
    auto m2 = text_to_rdmol(text, sample);
    BOOST_TEST(rdmol_to_text(*m2, sample) == text);

    // Check format auto-detect
    auto m3 = text_to_rdmol(text);
    BOOST_TEST(rdmol_to_text(*m3, sample) == text);

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol("garbage", sample),
                                    std::invalid_argument,
                                    "Failed to parse text");
}

BOOST_DATA_TEST_CASE(test_bypass_sanitization,
                     boost::unit_test::data::make(TEXT_FORMATS))
{
    if (sample == Format::SMARTS) {
        return; // skip SMARTS, which doesn't have sanitize options
    }
    // Create an unsanitized mol with a pentavalent C...
    auto mol = text_to_rdmol("C[C](C)(C)(C)C", Format::SMILES);
    auto text = rdmol_to_text(*mol, sample);

    // Make sure roundtripping preserves the poor chemistry
    auto m2 = text_to_rdmol(text, sample);
    BOOST_TEST(rdmol_to_text(*m2, sample) == text);
}

BOOST_AUTO_TEST_CASE(test_invalid_input)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol("garbage"),
                                    std::invalid_argument,
                                    "Unable to determine format");
}

BOOST_AUTO_TEST_CASE(test_roundtrip_smarts)
{
    // roundtrip through SMARTS
    std::string smarts = "cOc";
    auto mol = text_to_rdmol(smarts, Format::SMARTS);
    BOOST_TEST(rdmol_to_text(*mol, Format::SMARTS) == "cOc");

    // Read in as SMILES
    mol = text_to_rdmol(smarts, Format::SMILES);
    BOOST_TEST(rdmol_to_text(*mol, Format::SMARTS) == "[#6]-[#8]-[#6]");

    // Auto-detect reads as unsanitized SMILES
    mol = text_to_rdmol(smarts);
    BOOST_TEST(rdmol_to_text(*mol, Format::SMARTS) == "[#6]-[#8]-[#6]");
}

BOOST_AUTO_TEST_CASE(test_INCHI_KEY)
{
    auto mol = text_to_rdmol("C[C](C)(C)(C)C", Format::SMILES);
    auto text = rdmol_to_text(*mol, Format::INCHI_KEY);
    BOOST_TEST(text == "PJNAWCYHQGIPJJ-UHFFFAOYSA-N");

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol(text), std::invalid_argument,
                                    "Unable to determine format");

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol(text, Format::INCHI_KEY),
                                    std::invalid_argument,
                                    "Cannot read from INCHI_KEY");
}

BOOST_DATA_TEST_CASE(test_reactions_roundtrip,
                     boost::unit_test::data::make(REACTION_TEXT_FORMATS))
{
    std::string smiles = "CC(=O)O.OCC>>CC(=O)OCC";
    auto reaction = text_to_reaction(smiles, Format::SMILES);
    auto text = reaction_to_text(*reaction, sample);

    // Check roundtripping
    auto reaction2 = text_to_reaction(text, sample);
    BOOST_TEST(reaction_to_text(*reaction2, sample) == text);

    // Confirm reaction formats aren't picked up and fail
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol(smiles),
                                    std::invalid_argument,
                                    "Unable to determine format");

    // General failure
    std::string regular_smiles = "CC";
    auto reaction4 = text_to_reaction(smiles, Format::SMILES);
    reaction_to_text(*reaction4, sample);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_reaction(regular_smiles),
                                    std::invalid_argument,
                                    "Unable to determine format");
}

BOOST_DATA_TEST_CASE(testDoNotForceKekulizationOnExport,
                     boost::unit_test::data::make(TEXT_FORMATS))
{
    // SKETCH-1416

    // This SMILES cannot be kekulized: N atom's valence is ambiguous,
    // an explicit H is required to disambiguate.
    auto mol = text_to_rdmol("c1ccnc1", Format::SMILES);
    BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 5);

    if (sample == Format::INCHI || sample == Format::PDB) {
        // INCHI and PDB force kekulization, so we expect them to throw
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(
            rdmol_to_text(*mol, sample), RDKit::KekulizeException,
            "Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4");
    } else {
        // any other export format should not throw
        rdmol_to_text(*mol, sample);
    }
}

BOOST_DATA_TEST_CASE(testExportPreservesKekulizationState,
                     boost::unit_test::data::make(TEXT_FORMATS))
{
    // SKETCH-1416

    auto mol = text_to_rdmol("c1ccccc1", Format::SMILES);
    BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 6);

    // Aromatic mol
    {
        std::map<Format, std::string> references{
            {Format::SMILES, "c1ccccc1"},
            {Format::EXTENDED_SMILES, "c1ccccc1"},
            {Format::SMARTS, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"},
            {Format::MDL_MOLV2000, R"CTAB(
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
)CTAB"},
            {Format::MDL_MOLV3000, R"CTAB(
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 5 6
M  V30 6 4 6 1
)CTAB"},
            // INCHI is always kekulized
            {Format::INCHI, "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"},
            // PDB is always kekulized
            {Format::PDB,
             R"PDB(
CONECT    1    2    2    6
CONECT    2    3
CONECT    3    4    4
CONECT    4    5
CONECT    5    6    6
)PDB"},
        };

        auto molblock_out = rdmol_to_text(*mol, sample);
        BOOST_REQUIRE_NE(molblock_out.find(references[sample]),
                         std::string::npos);
    }

    // Kekulized mol
    RDKit::MolOps::Kekulize(*mol);
    {
        std::map<Format, std::string> references{
            {Format::SMILES, "C1=CC=CC=C1"},
            {Format::EXTENDED_SMILES, "C1=CC=CC=C1"},
            {Format::SMARTS, "[#6]1=[#6]-[#6]=[#6]-[#6]=[#6]-1"},
            {Format::MDL_MOLV2000, R"CTAB(
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
)CTAB"},
            {Format::MDL_MOLV3000, R"CTAB(
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
)CTAB"},
            // INCHI is always kekulized
            {Format::INCHI, "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"},
            // PDB is always kekulized
            {Format::PDB,
             R"PDB(
CONECT    1    2    2    6
CONECT    2    3
CONECT    3    4    4
CONECT    4    5
CONECT    5    6    6
)PDB"},
        };

        auto molblock_out = rdmol_to_text(*mol, sample);
        BOOST_REQUIRE_NE(molblock_out.find(references[sample]),
                         std::string::npos);
    }
}