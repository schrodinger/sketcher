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
