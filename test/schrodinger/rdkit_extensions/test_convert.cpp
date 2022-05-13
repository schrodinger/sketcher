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

    if (sample == Format::MDL_MOLV2000 || sample == Format::MDL_MOLV3000) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol("garbage", sample),
                                        std::invalid_argument,
                                        "Single molblock required; 0 present");
    } else {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(text_to_rdmol("garbage", sample),
                                        std::invalid_argument,
                                        "Failed to parse text");
    }
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

BOOST_AUTO_TEST_CASE(testAddEnhancedStereoToUngroupedMdlChiralAtoms1)
{
    // SKETCH-1453

    using RDKit::operator"" _ctab;

    for (auto chiral_flag : {1, 0}) {
        auto mol = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 Cl -3.179732 0.794933 0.000000 0
M  V30 2 C -1.742857 0.800000 0.000000 0
M  V30 3 H -1.020031 -0.441837 0.000000 0
M  V30 4 O -1.032939 2.039690 0.000000 0
M  V30 5 C -0.314287 0.797860 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END)MDL"_ctab;

        BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 4);

        mol->setProp(RDKit::common_properties::_MolFileChiralFlag, chiral_flag);

        add_enhanced_stereo_to_chiral_atoms(*mol);

        auto sgs = mol->getStereoGroups();
        BOOST_REQUIRE_EQUAL(sgs.size(), 1);

        auto sg = sgs.front();
        if (chiral_flag) {
            BOOST_CHECK(sg.getGroupType() ==
                        RDKit::StereoGroupType::STEREO_ABSOLUTE);
        } else {
            BOOST_CHECK(sg.getGroupType() ==
                        RDKit::StereoGroupType::STEREO_AND);
        }

        auto atoms = sg.getAtoms();
        BOOST_REQUIRE_EQUAL(atoms.size(), 1);
        // BOOST_CHECK_EQUAL(atoms.front(), mol->getAtomWithIdx(1));
    }
}

BOOST_AUTO_TEST_CASE(testAddEnhancedStereoToUngroupedMdlChiralAtoms2)
{
    // SKETCH-1453
    using RDKit::operator"" _ctab;

    struct test_ref {
        int chiral_flag;
        unsigned total_stereo_groups;
        unsigned abs_atoms_in_group;
        unsigned num_and_groups;
    };

    std::vector<test_ref> references = {
        test_ref{1, 2, 2, 1}, // ABS (2 atoms) + AND1 (1 atom)
        test_ref{0, 3, 1, 2}, // ABS (1 atom) + AND1 (1 atom) + AND2 (1 atom)
    };

    for (auto ref : references) {
        auto mol = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F -2.216000 -2.382857 0.000000 0
M  V30 2 C -2.218857 -0.954286 0.000000 0
M  V30 3 C -0.983143 -0.237429 0.000000 0
M  V30 4 C 0.255429 -0.949429 0.000000 0
M  V30 5 Cl 0.258286 -2.377714 0.000000 0
M  V30 6 C -3.457429 -0.242571 0.000000 0
M  V30 7 C -0.986000 1.191143 0.000000 0
M  V30 8 C 1.491143 -0.232571 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 2 6 CFG=1
M  V30 6 1 3 7 CFG=1
M  V30 7 1 4 8 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 MDLV30/STERAC1 ATOMS=(1 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)MDL"_ctab;

        BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 8);

        mol->setProp(RDKit::common_properties::_MolFileChiralFlag,
                     ref.chiral_flag);

        add_enhanced_stereo_to_chiral_atoms(*mol);

        unsigned abs_group_count{0};
        unsigned and_group_count{0};
        auto sgs = mol->getStereoGroups();

        BOOST_REQUIRE_EQUAL(sgs.size(), ref.total_stereo_groups);

        for (auto sg : sgs) {
            if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
                ++abs_group_count;
                BOOST_CHECK_EQUAL(sg.getAtoms().size(), ref.abs_atoms_in_group);
            } else if (sg.getGroupType() ==
                       RDKit::StereoGroupType::STEREO_AND) {
                ++and_group_count;
                BOOST_CHECK_EQUAL(sg.getAtoms().size(), 1);
            } else {
                BOOST_FAIL("unexpected Stereo Group Type");
            }
        }

        BOOST_CHECK_EQUAL(abs_group_count, 1);
        BOOST_CHECK_EQUAL(and_group_count, ref.num_and_groups);
    }
}
