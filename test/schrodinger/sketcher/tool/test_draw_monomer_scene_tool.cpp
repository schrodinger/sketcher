#define BOOST_TEST_MODULE test_draw_monomer_scene_tool

#include <memory>
#include <string>
#include <vector>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Helper that owns a monomer atom and its graphics item, and creates
 * UnboundMonomericAttachmentPointItem children for testing. The attachment
 * point items are Qt children of the monomer item, so they are deleted when the
 * monomer item is deleted.
 */
struct TestMonomer {
    std::unique_ptr<RDKit::Atom> atom;
    std::unique_ptr<AbstractMonomerItem> item;
    Fonts fonts;

    TestMonomer(const std::string& res_name,
                const rdkit_extensions::ChainType chain_type)
    {
        auto chain_id = rdkit_extensions::toString(chain_type) + "1";
        atom = rdkit_extensions::makeMonomer(res_name, chain_id, 1, false);
        item.reset(get_monomer_graphics_item(atom.get(), fonts));
    }

    /// Create an attachment point item with a numbered name (e.g. R1, R2)
    UnboundMonomericAttachmentPointItem* makeAP(int num)
    {
        UnboundAttachmentPoint ap{"R" + std::to_string(num), num, Direction::N};
        return new UnboundMonomericAttachmentPointItem(ap, item.get(), fonts);
    }

    /// Create an attachment point item with a custom name (num = -1)
    UnboundMonomericAttachmentPointItem* makeNamedAP(const std::string& name)
    {
        UnboundAttachmentPoint ap{name, ATTACHMENT_POINT_WITH_CUSTOM_NAME,
                                  Direction::N};
        return new UnboundMonomericAttachmentPointItem(ap, item.get(), fonts);
    }
};

/// Empty list returns nullptr
BOOST_AUTO_TEST_CASE(test_empty_list_returns_nullptr)
{
    std::vector<UnboundMonomericAttachmentPointItem*> empty;
    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::PEPTIDE, empty);
    BOOST_TEST(result == nullptr);
}

/// CHEM hovered: returns the attachment point with the lowest num
BOOST_AUTO_TEST_CASE(test_chem_hovered_returns_min_num)
{
    TestMonomer monomer("SMCC", rdkit_extensions::ChainType::CHEM);
    auto* ap3 = monomer.makeAP(3);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap3, ap1, ap2};

    auto* result = get_default_attachment_point(MonomerType::CHEM,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == ap1);
}

/// CHEM hovered: prefers numbered attachment points over custom-named ones
BOOST_AUTO_TEST_CASE(test_chem_hovered_prefers_numbered_over_custom_name)
{
    TestMonomer monomer("SMCC", rdkit_extensions::ChainType::CHEM);
    auto* ap_custom = monomer.makeNamedAP("pair");
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap_custom, ap2};

    auto* result = get_default_attachment_point(MonomerType::CHEM,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == ap2);

    // reverse the order of the attachment point list and make sure that we get
    // the same answer
    std::vector<UnboundMonomericAttachmentPointItem*> items_reversed{ap2,
                                                                     ap_custom};

    auto* result_reversed = get_default_attachment_point(
        MonomerType::CHEM, MonomerType::CHEM, items_reversed);
    BOOST_TEST(result_reversed == ap2);
}

/// CHEM hovered: returns custom-named attachment point when it's the only one
BOOST_AUTO_TEST_CASE(test_chem_hovered_returns_custom_name_when_only_option)
{
    TestMonomer monomer("SMCC", rdkit_extensions::ChainType::CHEM);
    auto* ap_custom = monomer.makeNamedAP("pair");
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap_custom};

    auto* result = get_default_attachment_point(MonomerType::CHEM,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == ap_custom);
}

/// PEPTIDE hovered + PEPTIDE tool: prefers R2 over R1 and R3
BOOST_AUTO_TEST_CASE(test_peptide_peptide_prefers_ap2)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    auto* ap3 = monomer.makeAP(3);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2, ap3};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == ap2);
}

/// PEPTIDE hovered + PEPTIDE tool: falls back to R1 when R2 is absent
BOOST_AUTO_TEST_CASE(test_peptide_peptide_falls_back_to_ap1)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap1 = monomer.makeAP(1);
    auto* ap3 = monomer.makeAP(3);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap3};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == ap1);
}

/// PEPTIDE hovered + PEPTIDE tool: falls back to R3 when R1 and R2 are absent
BOOST_AUTO_TEST_CASE(test_peptide_peptide_falls_back_to_ap3)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap3 = monomer.makeAP(3);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap3};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == ap3);
}

/// PEPTIDE hovered + PEPTIDE tool: returns nullptr when no preferred APs exist
BOOST_AUTO_TEST_CASE(test_peptide_peptide_returns_nullptr_when_no_preferred)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap4 = monomer.makeAP(4);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap4};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == nullptr);
}

/// PEPTIDE hovered + CHEM tool: returns R3
BOOST_AUTO_TEST_CASE(test_peptide_chem_returns_ap3)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    auto* ap3 = monomer.makeAP(3);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2, ap3};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == ap3);
}

/// PEPTIDE hovered + CHEM tool: returns nullptr when R3 is absent
BOOST_AUTO_TEST_CASE(test_peptide_chem_returns_nullptr_when_no_ap3)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == nullptr);
}

/// PEPTIDE hovered + unmatched tool type: returns nullptr
BOOST_AUTO_TEST_CASE(test_peptide_unmatched_tool_returns_nullptr)
{
    TestMonomer monomer("ALA", rdkit_extensions::ChainType::PEPTIDE);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2};

    auto* result = get_default_attachment_point(MonomerType::PEPTIDE,
                                                MonomerType::NA_BASE, items);
    BOOST_TEST(result == nullptr);
}

/// NA_BASE hovered + NA_BASE tool: returns the "pair" attachment point
BOOST_AUTO_TEST_CASE(test_na_base_na_base_returns_pair)
{
    TestMonomer monomer("A", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    auto* ap_pair = monomer.makeNamedAP("pair");
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap_pair};

    auto* result = get_default_attachment_point(MonomerType::NA_BASE,
                                                MonomerType::NA_BASE, items);
    BOOST_TEST(result == ap_pair);
}

/// NA_BASE hovered + NA_BASE tool: returns nullptr when no "pair" AP exists
BOOST_AUTO_TEST_CASE(test_na_base_na_base_returns_nullptr_when_no_pair)
{
    TestMonomer monomer("A", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(MonomerType::NA_BASE,
                                                MonomerType::NA_BASE, items);
    BOOST_TEST(result == nullptr);
}

/// NA_BASE hovered + CHEM tool: returns the "pair" attachment point
BOOST_AUTO_TEST_CASE(test_na_base_chem_returns_pair)
{
    TestMonomer monomer("A", rdkit_extensions::ChainType::RNA);
    auto* ap_pair = monomer.makeNamedAP("pair");
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap_pair};

    auto* result = get_default_attachment_point(MonomerType::NA_BASE,
                                                MonomerType::CHEM, items);
    BOOST_TEST(result == ap_pair);
}

/// NA_BASE hovered + NA_SUGAR tool: returns R1
BOOST_AUTO_TEST_CASE(test_na_base_na_sugar_returns_ap1)
{
    TestMonomer monomer("A", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    auto* ap_pair = monomer.makeNamedAP("pair");
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap_pair, ap1};

    auto* result = get_default_attachment_point(MonomerType::NA_BASE,
                                                MonomerType::NA_SUGAR, items);
    BOOST_TEST(result == ap1);
}

/// NA_BASE hovered + unmatched tool type: returns nullptr
BOOST_AUTO_TEST_CASE(test_na_base_unmatched_tool_returns_nullptr)
{
    TestMonomer monomer("A", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(
        MonomerType::NA_BASE, MonomerType::NA_PHOSPHATE, items);
    BOOST_TEST(result == nullptr);
}

/// NA_SUGAR hovered + NA_BASE tool: returns R3
BOOST_AUTO_TEST_CASE(test_na_sugar_na_base_returns_ap3)
{
    TestMonomer monomer("R", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    auto* ap3 = monomer.makeAP(3);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2, ap3};

    auto* result = get_default_attachment_point(MonomerType::NA_SUGAR,
                                                MonomerType::NA_BASE, items);
    BOOST_TEST(result == ap3);
}

/// NA_SUGAR hovered + NA_PHOSPHATE tool: prefers R2 over R1
BOOST_AUTO_TEST_CASE(test_na_sugar_na_phosphate_prefers_ap2)
{
    TestMonomer monomer("R", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2};

    auto* result = get_default_attachment_point(
        MonomerType::NA_SUGAR, MonomerType::NA_PHOSPHATE, items);
    BOOST_TEST(result == ap2);
}

/// NA_SUGAR hovered + NA_PHOSPHATE tool: falls back to R1 when R2 is absent
BOOST_AUTO_TEST_CASE(test_na_sugar_na_phosphate_falls_back_to_ap1)
{
    TestMonomer monomer("R", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(
        MonomerType::NA_SUGAR, MonomerType::NA_PHOSPHATE, items);
    BOOST_TEST(result == ap1);
}

/// NA_SUGAR hovered + unmatched tool type: returns nullptr
BOOST_AUTO_TEST_CASE(test_na_sugar_unmatched_tool_returns_nullptr)
{
    TestMonomer monomer("R", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(MonomerType::NA_SUGAR,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == nullptr);
}

/// NA_PHOSPHATE hovered + NA_SUGAR tool: prefers R2 over R1
BOOST_AUTO_TEST_CASE(test_na_phosphate_na_sugar_prefers_ap2)
{
    TestMonomer monomer("P", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    auto* ap2 = monomer.makeAP(2);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1, ap2};

    auto* result = get_default_attachment_point(MonomerType::NA_PHOSPHATE,
                                                MonomerType::NA_SUGAR, items);
    BOOST_TEST(result == ap2);
}

/// NA_PHOSPHATE hovered + NA_SUGAR tool: falls back to R1 when R2 is absent
BOOST_AUTO_TEST_CASE(test_na_phosphate_na_sugar_falls_back_to_ap1)
{
    TestMonomer monomer("P", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(MonomerType::NA_PHOSPHATE,
                                                MonomerType::NA_SUGAR, items);
    BOOST_TEST(result == ap1);
}

/// NA_PHOSPHATE hovered + unmatched tool type: returns nullptr
BOOST_AUTO_TEST_CASE(test_na_phosphate_unmatched_tool_returns_nullptr)
{
    TestMonomer monomer("P", rdkit_extensions::ChainType::RNA);
    auto* ap1 = monomer.makeAP(1);
    std::vector<UnboundMonomericAttachmentPointItem*> items{ap1};

    auto* result = get_default_attachment_point(MonomerType::NA_PHOSPHATE,
                                                MonomerType::PEPTIDE, items);
    BOOST_TEST(result == nullptr);
}

} // namespace sketcher
} // namespace schrodinger
