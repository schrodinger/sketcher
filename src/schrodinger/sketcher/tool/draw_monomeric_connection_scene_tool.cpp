#include "schrodinger/sketcher/tool/draw_monomeric_connection_scene_tool.h"

#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

DrawMonomericConnectionSceneTool::DrawMonomericConnectionSceneTool(
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    AbstractDrawMonomericConnectionSceneTool(fonts, scene, mol_model)
{
}

QPixmap DrawMonomericConnectionSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/monomeric_connection.svg");
}

std::optional<HintFragmentMonomerInfo>
DrawMonomericConnectionSceneTool::getHintFragmentMonomerInfoForDragEnd(
    const HintFragmentMonomerInfo& hint_start_monomer_info,
    const DragEndInfo& drag_end_info)
{
    if (std::holds_alternative<QPointF>(drag_end_info)) {
        auto end_pos = std::get<QPointF>(drag_end_info);
        // the MonomerType doesn't matter here (since we're not going to draw an
        // actual molecule) so we arbitrarily pick CHEM
        return HintFragmentMonomerInfo{{nullptr},
                                       MonomerType::CHEM,
                                       to_mol_xy(end_pos),
                                       "",
                                       NEW_MONOMER_FROM_INVALID_DRAG};
    } else {
        auto [hovered_monomer_item, drag_end_ap_item] =
            std::get<MonomerAndAPItems>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            hovered_monomer_item, drag_end_ap_item);
    }
}

void DrawMonomericConnectionSceneTool::updateHoveredUnboundAP(
    UnboundMonomericAttachmentPointItem* hovered_ap_item)
{
    // color the hovered AP black and color all of the others gray
    for (auto* ap_item : m_unbound_ap_items) {
        ap_item->setColor(ap_item == hovered_ap_item
                              ? m_unbound_ap_label_hover_color
                              : m_unbound_ap_label_color);
    }
}

UnboundMonomericAttachmentPointItem* DrawMonomericConnectionSceneTool::
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool /* no_default_if_click_should_mutate */) const
{
    auto [monomer, hovered_type] = getHoveredMonomerAndType();
    if (m_unbound_ap_items.empty()) {
        return nullptr;
    } else if (hovered_type == MonomerType::CHEM) {
        return find_min_attachment_point_by_num(m_unbound_ap_items);
    } else if (hovered_type == MonomerType::PEPTIDE) {
        return find_preferred_attachment_point_by_num(
            m_unbound_ap_items,
            {PeptideAP::C, PeptideAP::N, PeptideAP::SIDECHAIN});
    } else if (hovered_type == MonomerType::NA_BASE) {
        return find_attachment_point_with_name(m_unbound_ap_items,
                                               H_BOND_AP_MODEL_NAME);
    } else if (hovered_type == MonomerType::NA_SUGAR) {
        return find_preferred_attachment_point_by_num(
            m_unbound_ap_items,
            {NASugarAP::THREE_PRIME, NASugarAP::FIVE_PRIME});
    } else if (hovered_type == MonomerType::NA_PHOSPHATE) {
        return find_preferred_attachment_point_by_num(
            m_unbound_ap_items,
            {NAPhosphateAP::TO_NEXT_SUGAR, NAPhosphateAP::TO_PREV_SUGAR});
    }
    return nullptr;
}

} // namespace sketcher
} // namespace schrodinger
