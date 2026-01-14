#define BOOST_TEST_MODULE test_monomer_coloring

#include <string>

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

// Include the LocalMonomerDbFixture from rdkit_extensions tests
#include "../../rdkit_extensions/test_common.h"

using namespace schrodinger::sketcher;
using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_CASE(TestColorForNaturalMonomer)
{
    // Test that a natural amino acid monomer gets its expected color
    auto color = get_color_for_monomer("A", ChainType::PEPTIDE,
                                       AMINO_ACID_COLOR_BY_RES_NAME,
                                       DEFAULT_AA_BACKGROUND_COLOR);

    // "A" (Alanine) should get the ALIPHATIC color
    auto expected_color = MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC);
    BOOST_CHECK_EQUAL(color.name().toStdString(),
                      expected_color.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestColorForNonNaturalMonomerWithoutAnalog)
{
    // Test that a non-natural monomer without a natural analog gets the default
    // color
    auto color = get_color_for_monomer("DUMMY", ChainType::PEPTIDE,
                                       AMINO_ACID_COLOR_BY_RES_NAME,
                                       DEFAULT_AA_BACKGROUND_COLOR);

    // A non-existent monomer should get the default color
    BOOST_CHECK_EQUAL(color.name().toStdString(),
                      DEFAULT_AA_BACKGROUND_COLOR.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestColorForNonNaturalMonomerWithNaturalAnalog)
{
    // Test that a non-natural monomer with a natural analog gets the analog's
    // color

    // First, load a custom monomer that has "A" as its natural analog
    constexpr std::string_view custom_monomer_json =
        ("[{"
         "\"symbol\": \"CustomAla\","
         "\"polymer_type\": \"PEPTIDE\","
         "\"natural_analog\": \"A\"," // Natural analog is Alanine
         "\"smiles\": \"CCCC\","      // Dummy SMILES
         "\"core_smiles\": \"CCCC\","
         "\"name\": \"Custom Alanine Variant\","
         "\"monomer_type\": \"Backbone\","
         "\"author\": \"test\","
         "\"pdbcode\": \"CUST\""
         "}]");

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_monomer_json);

    // Verify the natural analog was set correctly
    auto custom_ala_analog =
        monomer_db.getNaturalAnalog("CustomAla", ChainType::PEPTIDE);
    BOOST_REQUIRE(custom_ala_analog.has_value());
    BOOST_CHECK_EQUAL(*custom_ala_analog, "A");

    // Now test that CustomAla gets the same color as its natural analog "A"
    auto color = get_color_for_monomer("CustomAla", ChainType::PEPTIDE,
                                       AMINO_ACID_COLOR_BY_RES_NAME,
                                       DEFAULT_AA_BACKGROUND_COLOR);

    // CustomAla should get the same color as "A" (ALIPHATIC)
    auto expected_color = MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC);
    BOOST_CHECK_EQUAL(color.name().toStdString(),
                      expected_color.name().toStdString());

    // Clean up
    monomer_db.resetMonomerDefinitions();
}

BOOST_AUTO_TEST_CASE(TestColorForNucleicAcidMonomerWithNaturalAnalog)
{
    // Test that a non-natural nucleic acid monomer with a natural analog gets
    // the analog's color

    // First, load a custom nucleic acid monomer that has "A" (Adenine) as its
    // natural analog
    constexpr std::string_view custom_monomer_json =
        ("[{"
         "\"symbol\": \"ModA\","
         "\"polymer_type\": \"RNA\","
         "\"natural_analog\": \"A\"," // Natural analog is Adenine
         "\"smiles\": \"NNNN\","      // Dummy SMILES
         "\"core_smiles\": \"NNNN\","
         "\"name\": \"Modified Adenine\","
         "\"monomer_type\": \"Backbone\","
         "\"author\": \"test\","
         "\"pdbcode\": \"MODA\""
         "}]");

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_monomer_json);

    // Verify the natural analog was set correctly
    auto mod_a_analog = monomer_db.getNaturalAnalog("ModA", ChainType::RNA);
    BOOST_REQUIRE(mod_a_analog.has_value());
    BOOST_CHECK_EQUAL(*mod_a_analog, "A");

    // Now test that ModA gets the same color as its natural analog "A"
    auto color = get_color_for_monomer(
        "ModA", ChainType::RNA, NUCLEIC_ACID_COLOR_BY_RES_NAME,
        MONOMER_COLOR_MAP.at(MonomerColorType::OTHER));

    // ModA should get the same color as "A" (ADENINE)
    auto expected_color = MONOMER_COLOR_MAP.at(MonomerColorType::ADENINE);
    BOOST_CHECK_EQUAL(color.name().toStdString(),
                      expected_color.name().toStdString());

    // Clean up
    monomer_db.resetMonomerDefinitions();
}
