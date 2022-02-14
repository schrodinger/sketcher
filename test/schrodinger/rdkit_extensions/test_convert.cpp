/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: text block <-> rdkit mol
 conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_convert

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/convert.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/sketcher.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::sketcher;
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
