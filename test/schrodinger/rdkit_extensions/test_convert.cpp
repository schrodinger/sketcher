/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: text block <-> rdkit mol
 conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE rdkit_extensions_convert

#include <map>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>

#include <fmt/format.h>

#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

namespace bdata = boost::unit_test::data;
using namespace schrodinger;
using namespace schrodinger::rdkit_extensions;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);
BOOST_TEST_DONT_PRINT_LOG_VALUE(RDKit::StereoGroupType);

BOOST_DATA_TEST_CASE(test_auto_detect,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    auto mol = to_rdkit("c1ccccc1", Format::SMILES);
    auto text = to_string(*mol, sample);

    // Check roundtripping and auto-detect
    for (auto tmp_mol : {to_rdkit(text, sample), to_rdkit(text)}) {
        BOOST_TEST(tmp_mol->getNumAtoms() == 6);
        BOOST_TEST(tmp_mol->getNumBonds() == 6);
        BOOST_TEST(to_string(*tmp_mol, sample) == text);
    }
}

BOOST_DATA_TEST_CASE(test_bypass_sanitization,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    // Create an unsanitized mol with a pentavalent C...
    auto mol = to_rdkit("C[C](C)(C)(C)C", Format::SMILES);

    if (sample == Format::XYZ) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(
            to_string(*mol, sample), RDKit::AtomValenceException,
            "Explicit valence for atom # 1 C, 5, is greater than permitted");
        return;
    }

    // Make sure roundtripping preserves the poor chemistry
    auto text = to_string(*mol, sample);
    auto m2 = to_rdkit(text, sample);
    BOOST_TEST(m2->getNumAtoms() == 6);
    BOOST_TEST(m2->getNumBonds() == 5);
    BOOST_TEST(to_string(*m2, sample) == text);
}

BOOST_DATA_TEST_CASE(test_parsing_error,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit("garbage", sample),
                                    std::invalid_argument,
                                    "Failed to parse text");
}

BOOST_AUTO_TEST_CASE(test_invalid_input)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit("garbage"), std::invalid_argument,
                                    "Unable to determine format");
}

BOOST_AUTO_TEST_CASE(test_roundtrip_smarts)
{
    // roundtrip through SMARTS
    std::string smarts = "cOc";
    auto mol = to_rdkit(smarts, Format::SMARTS);
    BOOST_TEST(to_string(*mol, Format::SMARTS) == smarts);

    // Read in as SMILES
    mol = to_rdkit(smarts, Format::SMILES);
    BOOST_TEST(to_string(*mol, Format::SMARTS) == "[#6]-[#8]-[#6]");

    // Auto-detect reads as unsanitized SMILES
    mol = to_rdkit(smarts);
    BOOST_TEST(to_string(*mol, Format::SMARTS) == "[#6]-[#8]-[#6]");

    // Test extended SMARTS
    smarts = "[#6]-[#8]-[#6]-* |$;;;_R1$|";
    mol = to_rdkit(smarts, Format::EXTENDED_SMARTS);
    BOOST_TEST(to_string(*mol, Format::EXTENDED_SMARTS) == smarts);
    // Auto-detect
    mol = to_rdkit(smarts);
    BOOST_TEST(to_string(*mol, Format::EXTENDED_SMARTS) == smarts);
}

BOOST_AUTO_TEST_CASE(test_INCHI_KEY)
{
    auto mol = to_rdkit("C[C](C)(C)(C)C", Format::SMILES);
    auto text = to_string(*mol, Format::INCHI_KEY);
    BOOST_TEST(text == "PJNAWCYHQGIPJJ-UHFFFAOYSA-N");

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit(text), std::invalid_argument,
                                    "Unable to determine format");

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit(text, Format::INCHI_KEY),
                                    std::invalid_argument,
                                    "Cannot read from INCHI_KEY");
}

BOOST_AUTO_TEST_CASE(test_MOL2)
{
    auto example = R"MOL(@<TRIPOS>MOLECULE
tmp.mol2: molecule 1
5   4    1
SMALL
NO_CHARGES


@<TRIPOS>ATOM
      1 C1         -0.1287    0.6148    0.0000 C.3       1 ****        0.0000
      2 H2         -0.7747    1.2608    0.6460 H         1 ****        0.0000
      3 H3         -0.7747   -0.0312   -0.6460 H         1 ****        0.0000
      4 H4          0.5173    1.2608   -0.6460 H         1 ****        0.0000
      5 H5          0.5173   -0.0312    0.6460 H         1 ****        0.0000
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
     3    1    4 1
     4    1    5 1
@<TRIPOS>SUBSTRUCTURE
     1 ****        1 GROUP             0       ****    0 ROOT)MOL";

    // Read both explicitly and via auto-detect
    for (auto mol : {to_rdkit(example, Format::MOL2), to_rdkit(example)}) {
        BOOST_TEST(mol->getNumAtoms() == 5);
        BOOST_TEST(mol->getNumBonds() == 4);
        BOOST_TEST(to_string(*mol, Format::SMILES) == "[H]C([H])([H])[H]");
    }

    // No write support
    auto mol = to_rdkit("C", Format::SMILES);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_string(*mol, Format::MOL2),
                                    std::invalid_argument,
                                    "Unsupported export format");
}

BOOST_AUTO_TEST_CASE(test_CDXML)
{
    auto example = R"CDXML(<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">
<CDXML BondLength="30.000000" LabelFont="3" CaptionFont="4">
    <fonttable id="1">
        <font id="2" charset="utf-8" name="Arial"/>
        <font id="3" charset="utf-8" name="Times New Roman"/>
    </fonttable>
    <colortable>
        <color r="1" g="1" b="1"/>
        <color r="0" g="0" b="0"/>
        <color r="1" g="0" b="0"/>
        <color r="1" g="1" b="0"/>
        <color r="0" g="1" b="0"/>
        <color r="0" g="1" b="1"/>
        <color r="0" g="0" b="1"/>
        <color r="1" g="0" b="1"/>
        <color r="0.5" g="0.5" b="0.5"/>
    </colortable>
    <page HeightPages="1" WidthPages="1">
        <fragment id="4">
            <n id="5" p="30.000000 30.000000">
                <t p="30.000000 30.000000" Justification="Center">
                    <s font="3" size="10" face="96">CH4</s>
                </t>
            </n>
        </fragment>
    </page>
</CDXML>)CDXML";

    // Read both explicitly and via auto-detect
    for (auto mol : {to_rdkit(example, Format::CDXML), to_rdkit(example)}) {
        BOOST_TEST(mol->getNumAtoms() == 1);
        BOOST_TEST(mol->getNumBonds() == 0);
        BOOST_TEST(to_string(*mol, Format::SMILES) == "C");
    }

    // No write support
    auto mol = to_rdkit("C", Format::SMILES);
    BOOST_REQUIRE_THROW(to_string(*mol, Format::CDXML), std::invalid_argument);
}

BOOST_DATA_TEST_CASE(test_reactions_roundtrip,
                     boost::unit_test::data::make(RXN_FORMATS))
{
    std::string smiles = "CC(=O)O.OCC>>CC(=O)OCC";
    auto reaction = to_rdkit_reaction(smiles, Format::SMILES);
    auto text = to_string(*reaction, sample);

    // Check roundtripping
    auto reaction2 = to_rdkit_reaction(text, sample);
    BOOST_TEST(to_string(*reaction2, sample) == text);

    // Confirm reaction formats aren't picked up and fail
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit(smiles), std::invalid_argument,
                                    "Unable to determine format");

    // General failure
    std::string regular_smiles = "CC";
    auto reaction4 = to_rdkit_reaction(smiles, Format::SMILES);
    to_string(*reaction4, sample);
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit_reaction(regular_smiles),
                                    std::invalid_argument,
                                    "Unable to determine format");
}

BOOST_DATA_TEST_CASE(testDoNotForceKekulizationOnExport,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    // SKETCH-1416

    // This SMILES cannot be kekulized: N atom's valence is ambiguous,
    // an explicit H is required to disambiguate.
    auto mol = to_rdkit("c1ccnc1", Format::SMILES);
    BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 5);

    // Formats that force kekulization are expected to throw
    if (sample == Format::INCHI || sample == Format::MAESTRO ||
        sample == Format::PDB || sample == Format::XYZ) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(
            to_string(*mol, sample), RDKit::KekulizeException,
            "Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4");
    } else {
        // any other export format should not throw
        to_string(*mol, sample);
    }
}

BOOST_DATA_TEST_CASE(testExportPreservesKekulizationState,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    // SKETCH-1416

    auto mol = to_rdkit("c1ccccc1", Format::SMILES);
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

        auto molblock_out = to_string(*mol, sample);
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

        auto molblock_out = to_string(*mol, sample);
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

        RDKit::translateChiralFlagToStereoGroups(*mol);

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

        RDKit::translateChiralFlagToStereoGroups(*mol);

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

BOOST_AUTO_TEST_CASE(testAddEnhancedStereoToUngroupedChiralAtomsSafety)
{
    // RDKit::translateChiralFlagToStereoGroups() should always add ABS groups
    // unless the input is coming from a mol block with chiral_flag=0.

    const auto smiles = "N[C@@H](CS)C(=O)O";
    {
        // non-to_rdkit inputs!
        const std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
        BOOST_REQUIRE_EQUAL(
            mol->hasProp(RDKit::common_properties::_MolFileChiralFlag), false);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        BOOST_REQUIRE_EQUAL(mol->getStereoGroups().size(), 0);

        mol->setProp(RDKit::common_properties::_MolFileChiralFlag, 1);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        const auto estg = mol->getStereoGroups();
        BOOST_REQUIRE_EQUAL(estg.size(), 1);
        BOOST_CHECK_EQUAL(estg.front().getGroupType(),
                          RDKit::StereoGroupType::STEREO_ABSOLUTE);
        BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), smiles);
    }
    {
        // to_rdkit SMILES input
        const auto mol = to_rdkit(smiles);
        BOOST_REQUIRE_EQUAL(
            mol->hasProp(RDKit::common_properties::_MolFileChiralFlag), false);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        BOOST_REQUIRE_EQUAL(mol->getStereoGroups().size(), 0);

        mol->setProp(RDKit::common_properties::_MolFileChiralFlag, 1);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        const auto estg = mol->getStereoGroups();
        BOOST_REQUIRE_EQUAL(estg.size(), 1);
        BOOST_CHECK_EQUAL(estg.front().getGroupType(),
                          RDKit::StereoGroupType::STEREO_ABSOLUTE);
        BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), smiles);
    }

    {
        const auto molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 1.299038 2.250000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 0.000000 0.000000 0.000000 0
M  V30 4 S -1.299038 0.750000 0.000000 0
M  V30 5 C 2.598076 -0.000000 0.000000 0
M  V30 6 O 2.598076 -1.500000 0.000000 0
M  V30 7 O 3.897114 0.750000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 2 5 6
M  V30 6 1 5 7
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
        const auto mol = to_rdkit(molblock);
        BOOST_REQUIRE_EQUAL(
            mol->hasProp(RDKit::common_properties::_MolFileChiralFlag), true);
        BOOST_CHECK_EQUAL(
            mol->getProp<int>(RDKit::common_properties::_MolFileChiralFlag), 0);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        const auto estg = mol->getStereoGroups();
        BOOST_REQUIRE_EQUAL(estg.size(), 1);
        BOOST_CHECK_EQUAL(estg.front().getGroupType(),
                          RDKit::StereoGroupType::STEREO_AND);

        // note that the CXSMILES is canonicalized (parity is reversed), but the
        // standard SMILES is not!
        const auto ps = RDKit::SmilesWriteParams();
        const auto flags =
            RDKit::SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS;
        BOOST_CHECK_EQUAL(RDKit::MolToCXSmiles(*mol, ps, flags),
                          "N[C@H](CS)C(=O)O |&1:1|");
        BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), "N[C@@H](CS)C(=O)O");
    }
    {
        const auto molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N 1.299038 2.250000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 0.000000 0.000000 0.000000 0
M  V30 4 S -1.299038 0.750000 0.000000 0
M  V30 5 C 2.598076 -0.000000 0.000000 0
M  V30 6 O 2.598076 -1.500000 0.000000 0
M  V30 7 O 3.897114 0.750000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 2 5 6
M  V30 6 1 5 7
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
        const auto mol = to_rdkit(molblock);
        BOOST_REQUIRE_EQUAL(
            mol->hasProp(RDKit::common_properties::_MolFileChiralFlag), true);
        BOOST_CHECK_EQUAL(
            mol->getProp<int>(RDKit::common_properties::_MolFileChiralFlag), 1);
        RDKit::translateChiralFlagToStereoGroups(*mol);
        const auto estg = mol->getStereoGroups();
        BOOST_REQUIRE_EQUAL(estg.size(), 1);
        BOOST_CHECK_EQUAL(estg.front().getGroupType(),
                          RDKit::StereoGroupType::STEREO_ABSOLUTE);
        BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), smiles);
    }
}

BOOST_AUTO_TEST_CASE(test_molTotValence_cleanup)
{
    // SHARED-8901

    auto molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -3.657143 -0.742857 0.000000 0 CHG=1 VAL=4
M  V30 END ATOM
M  V30 END CTAB
M  END)CTAB";

    auto mol = to_rdkit(molblock, Format::MDL_MOLV3000);
    BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 1);

    auto cxsmiles = to_string(*mol, Format::EXTENDED_SMILES);
    BOOST_TEST(cxsmiles.find("molTotValence") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(test_force_v2k)
{
    auto molblock = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.742857 1.028571 0.000000 0
M  V30 2 N -4.742857 2.457143 0.000000 0
M  V30 3 O -5.980036 0.314286 0.000000 0
M  V30 4 C -3.505678 0.314286 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  END)MDL";

    auto mol = to_rdkit(molblock);
    BOOST_TEST(mol->getNumAtoms() == 4);

    auto v3k = to_string(*mol, Format::MDL_MOLV3000);
    BOOST_TEST(v3k.find("V3000") != std::string::npos);
    auto v2k = to_string(*mol, Format::MDL_MOLV2000);
    BOOST_TEST(v2k.find("V2000") != std::string::npos);

    // V2K blows away all enhanced stereo!!!
    mol = to_rdkit(v2k);
    BOOST_TEST(mol->getStereoGroups().size() == 0);
}

BOOST_AUTO_TEST_CASE(test_attachment_point_coords)
{
    auto molblock = R"MDL(
  Mrv2216 10242223162D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.4167 2.9583 0 0 ATTCHPT=-1
M  V30 END ATOM
M  V30 END CTAB
M  END)MDL";

    auto mol = to_rdkit(molblock);
    auto conf = mol->getConformer();
    auto a0 = conf.getAtomPos(0);
    auto a1 = conf.getAtomPos(1);
    auto a2 = conf.getAtomPos(2);
    Eigen::Vector3d a(a0.x, a0.y, a0.z);
    Eigen::Vector3d b(a1.x, a1.y, a1.z);
    Eigen::Vector3d c(a2.x, a2.y, a2.z);
    // *C* is linear, the cross product is the zero vector
    BOOST_TEST((a - b).cross(a - c).norm() != 0.0);
}

BOOST_AUTO_TEST_CASE(test_atom_ring_queries)
{
    auto molblock = R""""(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.314286 0.828571 0.000000 0
M  V30 2 C -6.077107 1.542857 0.000000 0
M  V30 3 C -4.839927 0.828571 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) FIELDDISP="    0.0000    0.0000    DR    ALL  0 0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&R]"
M  V30 END SGROUP
M  V30 END CTAB
M  END)"""";

    auto mol = to_rdkit(molblock);
    auto roundtrip_mol = to_rdkit(to_string(*mol, Format::MDL_MOLV3000));

    // roundtripped mol should have the correct ring query atom
    auto at = roundtrip_mol->getAtomWithIdx(1);
    BOOST_REQUIRE(at->hasQuery());

    auto query_description = RDKit::describeQuery(at);
    BOOST_TEST(query_description.find("AtomInNRings") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(test_XYZ_import)
{
    std::string xyz_input = R"XYZ(7
CC(=O)[O-]; charge=-1; multiplicity=2
C     -1.239714   -0.709714    0.000000
C      0.000000   -0.000000    0.000000
O      0.005143    1.428571    0.000000
O      1.234571   -0.718571    0.000000
H     -1.060579   -1.693539    0.000000
H     -2.146426   -0.287964    0.000000
H     -2.031075   -1.321063    0.000000

)XYZ";

    // Check explicit read
    auto m1 = to_rdkit(xyz_input, Format::XYZ);
    BOOST_REQUIRE(m1 != nullptr);
    BOOST_TEST(m1->getNumAtoms() == 4);
    BOOST_TEST(m1->getNumBonds() == 3);
    // where the charge is parsed, but title/multiplicity info is lost
    BOOST_TEST(!m1->hasProp(RDKit::common_properties::_Name));
    BOOST_TEST(RDKit::MolOps::getFormalCharge(*m1) == -1);
    BOOST_TEST(!m1->hasProp("i_m_Spin_multiplicity"));

    // Check auto-detect
    auto m2 = to_rdkit(xyz_input, Format::XYZ);
    BOOST_REQUIRE(m2 != nullptr);

    // Test empty or invalid change lines
    auto empty_charge = R"XYZ(3
charge=
O      0.000000    0.000000    0.000000
H      1.000000    0.000000    0.000000
H     -0.333330   -0.942810    0.000000
)XYZ";
    auto invalid_charge = R"XYZ(3
charge=unknown
O      0.000000    0.000000    0.000000
H      1.000000    0.000000    0.000000
H     -0.333330   -0.942810    0.000000
)XYZ";
    for (auto bad_input : {empty_charge, invalid_charge}) {
        m1 = to_rdkit(bad_input, Format::XYZ);
        BOOST_REQUIRE(m1 != nullptr);
        BOOST_TEST(m1->getNumAtoms() == 1);
        BOOST_TEST(m1->getNumBonds() == 0);
        BOOST_TEST(RDKit::MolOps::getFormalCharge(*m1) == 0);
    }
}

BOOST_AUTO_TEST_CASE(test_XYZ_export)
{
    auto mol = to_rdkit("C[NH3+]");
    auto reference = R"XYZ(8
charge=1
C      0.706212    0.016657   -0.070197
N     -0.730887   -0.055064    0.057318
H      1.234241   -0.189593    0.886397
H      0.938111    1.079963   -0.325607
H      1.124059   -0.656270   -0.833721
H     -0.937185    0.111382    1.078681
H     -1.211908    0.671553   -0.537460
H     -1.122643   -0.978629   -0.255412
)XYZ";
    BOOST_TEST(to_string(*mol, Format::XYZ) == reference);
    mol->setProp(RDKit::common_properties::_Name, "C[NH3+]");
    reference = R"XYZ(8
C[NH3+]; charge=1
C      0.706212    0.016657   -0.070197
N     -0.730887   -0.055064    0.057318
H      1.234241   -0.189593    0.886397
H      0.938111    1.079963   -0.325607
H      1.124059   -0.656270   -0.833721
H     -0.937185    0.111382    1.078681
H     -1.211908    0.671553   -0.537460
H     -1.122643   -0.978629   -0.255412
)XYZ";
    BOOST_TEST(to_string(*mol, Format::XYZ) == reference);

    mol = to_rdkit("CC(=O)[O-]");
    reference = R"XYZ(7
charge=-1
C     -0.671981   -0.096609    0.012890
C      0.797773    0.013000   -0.239093
O      1.562004   -0.962883    0.028299
O      1.352785    1.268077    0.036242
H     -1.087340   -1.074207   -0.268495
H     -0.825770    0.169879    1.078053
H     -1.127471    0.682743   -0.647895
)XYZ";
    BOOST_TEST(to_string(*mol, Format::XYZ) == reference);
    mol->setProp(RDKit::common_properties::_Name, "CC(=O)[O-]");
    reference = R"XYZ(7
CC(=O)[O-]; charge=-1
C     -0.671981   -0.096609    0.012890
C      0.797773    0.013000   -0.239093
O      1.562004   -0.962883    0.028299
O      1.352785    1.268077    0.036242
H     -1.087340   -1.074207   -0.268495
H     -0.825770    0.169879    1.078053
H     -1.127471    0.682743   -0.647895
)XYZ";
    BOOST_TEST(to_string(*mol, Format::XYZ) == reference);

    mol = to_rdkit("O");
    mol->setProp("i_m_Spin_multiplicity", 2);
    reference = R"XYZ(3
multiplicity=2
O     -0.009295    0.396241   -0.000000
H     -0.786760   -0.201660   -0.000000
H      0.796055   -0.194581    0.000000
)XYZ";
    BOOST_TEST(to_string(*mol, Format::XYZ) == reference);
}

BOOST_AUTO_TEST_CASE(TestConvertingHELM)
{
    const std::string text{"PEPTIDE1{A.K.L}$$$$V2.0"};
    const auto mol = to_rdkit(text, Format::HELM);
    BOOST_TEST(to_string(*mol, Format::HELM) == text);
}

BOOST_DATA_TEST_CASE(test_cannot_be_smiles,
                     boost::unit_test::data::make({"[#6]-[#6]-[#7]-[#6]", "C~C",
                                                   "[13#6]-[#6]"}),
                     pseudo_smiles)
{
    // SHARED-10397, SKETCH-2190
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(
        to_rdkit(pseudo_smiles, Format::SMILES), std::invalid_argument,
        fmt::format("{} is not a valid SMILES", pseudo_smiles));

    auto mol = to_rdkit(pseudo_smiles, Format::AUTO_DETECT);
    auto smarts_mol = to_rdkit(pseudo_smiles, Format::SMARTS);
    BOOST_TEST(mol != nullptr);
    BOOST_TEST(smarts_mol != nullptr);
    BOOST_TEST(to_string(*mol, Format::SMARTS) ==
               to_string(*smarts_mol, Format::SMARTS));
}

BOOST_AUTO_TEST_CASE(test_smarts_no_radicals)
{
    // SKETCH-2190
    const auto smarts = "[CH3]~[CH2]~*";

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(
        to_rdkit(smarts, Format::SMILES), std::invalid_argument,
        fmt::format("{} is not a valid SMILES", smarts));

    auto mol = to_rdkit(smarts, Format::AUTO_DETECT);
    auto smarts_mol = to_rdkit(smarts, Format::SMARTS);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_REQUIRE(smarts_mol != nullptr);
    BOOST_TEST(to_string(*mol, Format::SMARTS) ==
               to_string(*smarts_mol, Format::SMARTS));

    for (auto atom : mol->atoms()) {
        BOOST_TEST(atom->getNumRadicalElectrons() == 0);
    }
}

/**
 * Make sure that we can parse extended SMILES and reaction SMILES/SMARTS
 * strings that use atom labels to indicate R-groups (see SHARED-10951).  Also
 * make sure that R-groups are exported to EXTENDED_SMILES using atom labels
 */
BOOST_AUTO_TEST_CASE(test_atom_mapping_rgroups)
{
    auto mol = to_rdkit("C* |$;_R1$|", Format::EXTENDED_SMILES);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_TEST(mol->getNumAtoms() == 2);
    auto r_group_num = get_r_group_number(mol->getAtomWithIdx(1));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 1);
    BOOST_TEST(to_string(*mol, Format::EXTENDED_SMILES) == "*C |$_R1;$|");

    mol = to_rdkit("C* |$;_R5$|", Format::EXTENDED_SMILES);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_TEST(mol->getNumAtoms() == 2);
    r_group_num = get_r_group_number(mol->getAtomWithIdx(1));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 5);
    BOOST_TEST(to_string(*mol, Format::EXTENDED_SMILES) == "*C |$_R5;$|");

    // test a simple reaction
    auto rxn = to_rdkit_reaction("C*.C*>>*C* |$;_R10;;_R20;_R10;;_R20$|");
    BOOST_REQUIRE(rxn != nullptr);
    auto reactant_1 = rxn->getReactants()[0];
    r_group_num = get_r_group_number(reactant_1->getAtomWithIdx(1));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 10);
    auto reactant_2 = rxn->getReactants()[1];
    r_group_num = get_r_group_number(reactant_2->getAtomWithIdx(1));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 20);
    auto product = rxn->getProducts()[0];
    r_group_num = get_r_group_number(product->getAtomWithIdx(0));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 10);
    r_group_num = get_r_group_number(product->getAtomWithIdx(2));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 20);

    // test the reaction from SHARED-10951
    rxn = to_rdkit_reaction(
        "*-C#N.Br>>*-N1-C(-*)=N-N=N-1 |$_R1;;;;_R2;;;_R1;;;;$|");
    BOOST_REQUIRE(rxn != nullptr);
    reactant_1 = rxn->getReactants()[0];
    r_group_num = get_r_group_number(reactant_1->getAtomWithIdx(0));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 1);
    product = rxn->getProducts()[0];
    r_group_num = get_r_group_number(product->getAtomWithIdx(0));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 2);
    r_group_num = get_r_group_number(product->getAtomWithIdx(3));
    BOOST_TEST(r_group_num.has_value());
    BOOST_TEST(*r_group_num == 1);
}

// Since we can convert monomeric mols to atomistic mols, we shouldn't
// have any issues converting them to all supported mol formats.
BOOST_DATA_TEST_CASE(test_converting_biologics_to_non_reaction_formats,
                     bdata::make(MOL_FORMATS), mol_format)
{
    using test_info_t = std::pair<std::string, Format>;
    for (auto [input, input_format] : std::array<test_info_t, 4>{{
             {">\nACG\n", Format::FASTA_DNA},
             {">\nACG\n", Format::FASTA_PEPTIDE},
             {">\nACG\n", Format::FASTA_RNA},
             {"PEPTIDE1{A.K.L}$$$$V2.0", Format::HELM},

         }}) {
        auto input_mol = to_rdkit(input, input_format);
        BOOST_CHECK_NO_THROW(to_string(*input_mol, mol_format));
    }
}

// Since we can convert atomistic mols to monomeric mols, we shouldn't
// have any issues converting them to sequence formats.
//
// NOTE: Choosing MAESTRO and PDB as the atomistic formats because they include
// residue information, which is required for constructing monomeric mols.
BOOST_DATA_TEST_CASE(test_converting_atomistic_biologics_to_seq_formats,
                     bdata::make(std::vector<Format>{
                         Format::MAESTRO,
                         Format::PDB,
                     }),
                     mol_format)
{
    auto monomeric_mol = to_rdkit("PEPTIDE1{A.C.G}$$$$V2.0", Format::HELM);
    auto atomistic_mol =
        to_rdkit(to_string(*monomeric_mol, mol_format), mol_format);
    BOOST_TEST(to_string(*atomistic_mol, Format::HELM) ==
               "PEPTIDE1{A.C.G}$$$$V2.0");
    BOOST_TEST(to_string(*atomistic_mol, Format::FASTA) == ">\nACG\n");
}

BOOST_DATA_TEST_CASE(test_converting_structures_with_wiggly_bonds,
                     bdata::make(std::vector<Format>{
                         Format::MDL_MOLV3000,
                         Format::EXTENDED_SMILES,
                     }),
                     test_format)
{
    using RDKit::common_properties::_MolFileBondCfg;

    // read structure with wiggly bonds
    auto mol = to_rdkit("CCC(C)N |w:2.3|", Format::EXTENDED_SMILES);

    // roundtrip mol through test format
    auto roundtripped_mol = to_rdkit(to_string(*mol, test_format), test_format);

    auto bond_config = 0u;
    // this should be the bond between the N and non-branch C
    auto test_bond = roundtripped_mol->getBondWithIdx(3);
    BOOST_REQUIRE(test_bond->getPropIfPresent(_MolFileBondCfg, bond_config));
    // Wiggly bonds should have config == 2
    BOOST_TEST(bond_config == 2);
    // Wiggly bonds should have UNKNOWN Direction
    BOOST_TEST(test_bond->getBondDir() == RDKit::Bond::UNKNOWN);
}

BOOST_DATA_TEST_CASE(TestCombiningAbsoluteEnhancedStereoGroups,
                     (bdata::make(std::vector<std::string>{
                          "STEABS",
                      }) ^
                      bdata::make(std::vector<std::string>{
                          "|a:2,4|",
                      })) *
                         (bdata::make(std::vector<std::string>{
                              "3",
                              "4",
                          }) ^
                          bdata::make(std::vector<std::string>{
                              "4",
                              "3",
                          })),
                     test_stereo_group_type, expected_stereo_group_extension,
                     test_atom1, test_atom2)
{
    auto molblock = fmt::format(R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.092597 -1.392509 0.000000 0
M  V30 2 C -2.683821 -0.875482 0.000000 0
M  V30 3 C -4.347972 -2.872425 0.000000 0
M  V30 4 C -3.194885 -3.830086 0.000000 0
M  V30 5 C -1.787157 -3.315618 0.000000 0
M  V30 6 C -3.451094 -5.309555 0.000000 0
M  V30 7 C -4.859281 -5.826684 0.000000 0
M  V30 8 C -5.755263 -3.391590 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 4 3
M  V30 4 1 5 4
M  V30 5 1 4 6 CFG=1
M  V30 6 1 6 7
M  V30 7 1 3 8 CFG=3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/{0} ATOMS=(1 {1})
M  V30 MDLV30/{0} ATOMS=(1 {2})
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$)MDL",
                                test_stereo_group_type, test_atom1, test_atom2);
    auto mol = to_rdkit(molblock);

    // we should only have one stereo group
    BOOST_TEST(mol->getStereoGroups().size() == 1);

    // clear coordinates to make the cxsmiles output cleaner
    mol->clearConformers();
    auto cxsmiles = to_string(*mol, Format::EXTENDED_SMILES);

    BOOST_TEST(cxsmiles.ends_with(expected_stereo_group_extension));
}

BOOST_AUTO_TEST_CASE(TestCombineStereoGroupsVectorSizingIssue)
{
    // SHARED-11993: Previous workaround for combining stereogroups lead to
    // a dangling pointer, ensure that ABS groups are correctly combined
    {
        auto molblock = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.321088 1.959235 0.000000 0
M  V30 2 C 1.087688 2.476262 0.000000 0
M  V30 3 C -0.576463 0.479319 0.000000 0
M  V30 4 C 0.576624 -0.478342 0.000000 0
M  V30 5 C 1.984352 0.036126 0.000000 0
M  V30 6 C 0.320415 -1.957811 0.000000 0
M  V30 7 C -1.087772 -2.474940 0.000000 0
M  V30 8 C -1.983754 -0.039846 0.000000 0
M  V30 9 C 2.239331 1.515160 0.000000 0
M  V30 10 N 1.344205 3.954165 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 4 3
M  V30 4 1 5 4
M  V30 5 1 4 6 CFG=1
M  V30 6 1 6 7
M  V30 7 1 3 8 CFG=3
M  V30 8 1 2 9
M  V30 9 1 2 10 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 3)
M  V30 MDLV30/STERAC1 ATOMS=(1 4)
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$)MDL";
        auto mol = to_rdkit(molblock);

        // we should only have two stereo groups now == 1 ABS, 1 racemic
        BOOST_TEST(mol->getStereoGroups().size() == 2);

        // clear coordinates to make the cxsmiles output cleaner
        mol->clearConformers();
        auto cxsmiles = to_string(*mol, Format::EXTENDED_SMILES);
        BOOST_TEST(cxsmiles == "CC[C@H](C)[C@H](C)C[C@H](C)N |a:4,7,&1:2|");
    }

    {
        auto molblock = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.321088 1.959235 0.000000 0
M  V30 2 C 1.087688 2.476262 0.000000 0
M  V30 3 C -0.576463 0.479319 0.000000 0
M  V30 4 C 0.576624 -0.478342 0.000000 0
M  V30 5 C 1.984352 0.036126 0.000000 0
M  V30 6 C 0.320415 -1.957811 0.000000 0
M  V30 7 C -1.087772 -2.474940 0.000000 0
M  V30 8 C -1.983754 -0.039846 0.000000 0
M  V30 9 C 2.239331 1.515160 0.000000 0
M  V30 10 N 1.344205 3.954165 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 4 3
M  V30 4 1 5 4
M  V30 5 1 4 6 CFG=1
M  V30 6 1 6 7
M  V30 7 1 3 8 CFG=3
M  V30 8 1 2 9
M  V30 9 1 2 10 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 3)
M  V30 MDLV30/STEREL1 ATOMS=(1 4)
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$)MDL";
        auto mol = to_rdkit(molblock);

        // we should only have two stereo groups now == 1 ABS, 1 relative
        BOOST_TEST(mol->getStereoGroups().size() == 2);

        // clear coordinates to make the cxsmiles output cleaner
        mol->clearConformers();
        auto cxsmiles = to_string(*mol, Format::EXTENDED_SMILES);
        BOOST_TEST(cxsmiles == "CC[C@H](C)[C@H](C)C[C@H](C)N |a:4,7,o1:2|");
    }
}