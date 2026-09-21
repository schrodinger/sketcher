#define BOOST_TEST_MODULE test_amino_acid_item

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/amino_acid_item.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

#include "../qapplication_required_fixture.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestAminoAcidItem : public AminoAcidItem
{
  public:
    using AminoAcidItem::AminoAcidItem;
    using AbstractMonomerItem::m_main_label_is_elided;
    using AbstractMonomerItem::m_main_label_paint_text;
    using AbstractMonomerItem::m_main_label_text;
};

BOOST_AUTO_TEST_CASE(test_long_label_is_shortened_for_fading)
{
    Fonts fonts;
    AtomDisplaySettings atom_display_settings;
    BondDisplaySettings bond_display_settings;

    auto monomer =
        rdkit_extensions::makeMonomer("FADINGX", "PEPTIDE1", 1, false);
    TestAminoAcidItem item(monomer.get(), fonts, atom_display_settings,
                           bond_display_settings);

    BOOST_TEST(item.m_main_label_text.toStdString() == "FADIN");
    BOOST_TEST(item.m_main_label_paint_text.toStdString() == "FADING");
    BOOST_TEST(item.m_main_label_is_elided);
    BOOST_TEST(elide_text("FADIN").toStdString() == "FADIN");

    auto six_character_monomer =
        rdkit_extensions::makeMonomer("FADING", "PEPTIDE1", 1, false);
    TestAminoAcidItem six_character_item(
        six_character_monomer.get(), fonts, atom_display_settings,
        bond_display_settings);
    BOOST_TEST(six_character_item.m_main_label_text.toStdString() == "FADING");
    BOOST_TEST(six_character_item.m_main_label_paint_text.toStdString() ==
               "FADING");
    BOOST_TEST(!six_character_item.m_main_label_is_elided);
}

/**
 * Make sure that the tool tip for a standard amino acid is empty, while the
 * tool tip for a SMILES monomer contains an image
 */
BOOST_AUTO_TEST_CASE(test_tool_tips)
{
    Fonts fonts;
    AtomDisplaySettings atom_display_settings;
    BondDisplaySettings bond_display_settings;

    auto alanine = rdkit_extensions::makeMonomer("A", "PEPTIDE1", 1, false);
    AminoAcidItem alanine_item(alanine.get(), fonts, atom_display_settings,
                               bond_display_settings);
    BOOST_TEST(alanine_item.toolTip().isEmpty());

    auto smiles_monomer =
        rdkit_extensions::makeMonomer("CC", "PEPTIDE1", 2, true);
    AminoAcidItem smiles_monomer_item(smiles_monomer.get(), fonts,
                                      atom_display_settings,
                                      bond_display_settings);
    BOOST_TEST(smiles_monomer_item.toolTip().startsWith("<img"));
}

} // namespace sketcher
} // namespace schrodinger
