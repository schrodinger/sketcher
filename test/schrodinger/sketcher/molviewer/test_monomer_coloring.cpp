#define BOOST_TEST_MODULE test_monomer_coloring

#include <string>

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/chem_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/nucleic_acid_base_item.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"

// Include the LocalMonomerDbFixture from rdkit_extensions tests
#include "../../rdkit_extensions/test_common.h"
#include "../qapplication_required_fixture.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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

// Test fixture that creates monomer graphics items from a HELM string and
// keeps the mol and fonts alive for the duration of the test.
struct HelmItemsFixture {
    AtomDisplaySettings atom_display_settings;
    BondDisplaySettings bond_display_settings;
    Fonts fonts;
    boost::shared_ptr<RDKit::RWMol> mol;
    std::vector<QGraphicsItem*> all_items;
    std::unordered_map<const RDKit::Atom*, QGraphicsItem*> atom_items;
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*> bond_items;
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*> secondary_items;
    std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*> sgroup_items;

    HelmItemsFixture(const std::string& helm)
    {
        mol = to_rdkit(helm);
        prepare_mol(*mol);
        for (auto* atom : mol->atoms()) {
            set_atom_monomeric(atom);
        }
        std::tie(all_items, atom_items, bond_items, secondary_items,
                 sgroup_items) =
            create_graphics_items_for_mol(
                mol.get(), fonts, atom_display_settings, bond_display_settings);
    }

    ~HelmItemsFixture()
    {
        for (auto* item : all_items) {
            delete item;
        }
    }

    QGraphicsItem* atom_item(unsigned int idx)
    {
        return atom_items.at(mol->getAtomWithIdx(idx));
    }

    QGraphicsItem* bond_item(unsigned int idx)
    {
        return bond_items.at(mol->getBondWithIdx(idx));
    }
};

BOOST_AUTO_TEST_CASE(TestDarkModeStandardAminoAcidBorderColors)
{
    // Atom 0 = "A" (Alanine), a standard amino acid
    HelmItemsFixture f("PEPTIDE1{A.G}$$$$V2.0");

    auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(f.atom_item(0));
    BOOST_REQUIRE(monomer_item != nullptr);

    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      STANDARD_AA_BORDER_COLOR.name().toStdString());

    monomer_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      STANDARD_AA_BORDER_COLOR_DARK_BG.name().toStdString());

    monomer_item->setDarkMode(false);
    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      STANDARD_AA_BORDER_COLOR.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeDAAminoAcidBorderColors)
{
    // "dA" starts with 'd', so it's classified as a D-amino acid
    HelmItemsFixture f("PEPTIDE1{[dA].A}$$$$V2.0");

    auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(f.atom_item(0));
    BOOST_REQUIRE(monomer_item != nullptr);

    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      D_AA_BORDER_COLOR.name().toStdString());

    monomer_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      D_AA_BORDER_COLOR_DARK_BG.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeOtherAminoAcidBorderColors)
{
    // "meA" is not in the standard map and doesn't start with 'd', so it's
    // OTHER
    HelmItemsFixture f("PEPTIDE1{[meA].A}$$$$V2.0");

    auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(f.atom_item(0));
    BOOST_REQUIRE(monomer_item != nullptr);

    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      OTHER_AA_BORDER_COLOR.name().toStdString());

    monomer_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      OTHER_AA_BORDER_COLOR_DARK_BG.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeNucleicAcidBaseBorderColors)
{
    // RNA1{R(A)P}: atom 0 = sugar (R), atom 1 = base (A), atom 2 = phosphate
    HelmItemsFixture f("RNA1{R(A)P.R(C)P.R(G)P}$$$$V2.0");

    auto* base_item = dynamic_cast<NucleicAcidBaseItem*>(f.atom_item(1));
    BOOST_REQUIRE(base_item != nullptr);

    BOOST_CHECK_EQUAL(base_item->getBorderColor().name().toStdString(),
                      STANDARD_NA_BORDER_COLOR.name().toStdString());

    base_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(base_item->getBorderColor().name().toStdString(),
                      STANDARD_NA_BORDER_COLOR_DARK_BG.name().toStdString());

    base_item->setDarkMode(false);
    BOOST_CHECK_EQUAL(base_item->getBorderColor().name().toStdString(),
                      STANDARD_NA_BORDER_COLOR.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeChemMonomerBorderColors)
{
    HelmItemsFixture f("CHEM1{*}$$$$V2.0");

    auto* chem_item = dynamic_cast<ChemMonomerItem*>(f.atom_item(0));
    BOOST_REQUIRE(chem_item != nullptr);

    BOOST_CHECK_EQUAL(chem_item->getBorderColor().name().toStdString(),
                      CHEM_MONOMER_BORDER_COLOR.name().toStdString());

    chem_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(chem_item->getBorderColor().name().toStdString(),
                      CHEM_MONOMER_BORDER_COLOR_DARK_BG.name().toStdString());

    chem_item->setDarkMode(false);
    BOOST_CHECK_EQUAL(chem_item->getBorderColor().name().toStdString(),
                      CHEM_MONOMER_BORDER_COLOR.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeConnectorColors)
{
    // Bond 0 connects atom 0 (A) to atom 1 (G)
    HelmItemsFixture f("PEPTIDE1{A.G}$$$$V2.0");

    auto* connector_item = dynamic_cast<MonomerConnectorItem*>(f.bond_item(0));
    BOOST_REQUIRE(connector_item != nullptr);

    BOOST_CHECK_EQUAL(connector_item->getConnectorColor().name().toStdString(),
                      AA_LINEAR_CONNECTOR_COLOR.name().toStdString());

    connector_item->setDarkMode(true);
    BOOST_CHECK_EQUAL(connector_item->getConnectorColor().name().toStdString(),
                      AA_LINEAR_CONNECTOR_COLOR_DARK_BG.name().toStdString());

    connector_item->setDarkMode(false);
    BOOST_CHECK_EQUAL(connector_item->getConnectorColor().name().toStdString(),
                      AA_LINEAR_CONNECTOR_COLOR.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeMonomerPersistsAfterUpdateCachedData)
{
    // Verify that dark mode border colors survive updateCachedData() (e.g. on
    // zoom).
    HelmItemsFixture f("PEPTIDE1{A.G}$$$$V2.0");

    auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(f.atom_item(0));
    BOOST_REQUIRE(monomer_item != nullptr);

    monomer_item->setDarkMode(true);
    monomer_item->updateCachedData();

    BOOST_CHECK_EQUAL(monomer_item->getBorderColor().name().toStdString(),
                      STANDARD_AA_BORDER_COLOR_DARK_BG.name().toStdString());
}

BOOST_AUTO_TEST_CASE(TestDarkModeConnectorPersistsAfterUpdateCachedData)
{
    // Verify that dark mode connector colors survive updateCachedData()
    HelmItemsFixture f("PEPTIDE1{A.G}$$$$V2.0");

    auto* connector_item = dynamic_cast<MonomerConnectorItem*>(f.bond_item(0));
    BOOST_REQUIRE(connector_item != nullptr);

    connector_item->setDarkMode(true);
    connector_item->updateCachedData();

    BOOST_CHECK_EQUAL(connector_item->getConnectorColor().name().toStdString(),
                      AA_LINEAR_CONNECTOR_COLOR_DARK_BG.name().toStdString());
}
