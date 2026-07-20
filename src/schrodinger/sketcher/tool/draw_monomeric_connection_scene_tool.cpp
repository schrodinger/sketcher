#include "schrodinger/sketcher/tool/draw_monomeric_connection_scene_tool.h"

#include <optional>
#include <utility>
#include <variant>

#include <QGraphicsLineItem>

#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

DrawMonomericConnectionSceneTool::DrawMonomericConnectionSceneTool(
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    // the residue name and chain type are irrelevant since this tool will never
    // create a new residue, so just pass dummy values
    AbstractDrawMonomerOrMonomericConnectionSceneTool(
        "", rdkit_extensions::ChainType::CHEM, fonts, scene, mol_model)
{
}

QPixmap DrawMonomericConnectionSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/monomeric_connection.svg");
}

void DrawMonomericConnectionSceneTool::updateColorsAfterBackgroundColorChange(
    bool is_dark_mode)
{
    AbstractDrawMonomerOrMonomericConnectionSceneTool::
        updateColorsAfterBackgroundColorChange(is_dark_mode);
    m_invalid_drag_pen.setColor(is_dark_mode
                                    ? CONNECTION_TOOL_INVALID_DRAG_COLOR_DARK_BG
                                    : CONNECTION_TOOL_INVALID_DRAG_COLOR);
    m_unbound_ap_label_hover_color =
        is_dark_mode ? CONNECTION_TOOL_AP_LABEL_HOVER_COLOR_DARK_BG
                     : CONNECTION_TOOL_AP_LABEL_HOVER_COLOR;
}

bool DrawMonomericConnectionSceneTool::clickShouldMutate(
    const RDKit::Atom* /* monomer */,
    const MonomerType /* monomer_type */) const
{
    return false;
}

std::optional<HintFragmentMonomerInfo>
DrawMonomericConnectionSceneTool::getHintFragmentMonomerInfoForDragStart()
{
    if (m_drag_start_monomer_item == nullptr) {
        // the drag started over empty space, so we do nothing (unlike
        // DrawMonomerSceneTool, which would add a new monomer here)
        return std::nullopt;
    } else if (!m_drag_start_ap_model_name.empty()) {
        // the drag started over a monomer and it has an available attachment
        // point
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            m_drag_start_monomer_item, m_drag_start_ap_model_name);
    } else {
        // the drag started over a monomer, but that monomer has no available
        // unbound attachment points so we can't drag from it
        return std::nullopt;
    }
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

DragEndInfo DrawMonomericConnectionSceneTool::getDragEndInfoForNonMonomerPos(
    const QPointF& cur_scene_pos) const
{
    return cur_scene_pos;
}

void DrawMonomericConnectionSceneTool::clearHintFragmentItem()
{
    AbstractDrawMonomerOrMonomericConnectionSceneTool::clearHintFragmentItem();
    delete m_invalid_drag_item;
    m_invalid_drag_item = nullptr;
}

void DrawMonomericConnectionSceneTool::createHintFragmentItem(
    HintFragmentMonomerInfo& monomer_one_info,
    HintFragmentMonomerInfo& monomer_two_info)
{
    if (monomer_one_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG &&
        monomer_two_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG) {
        AbstractDrawMonomerOrMonomericConnectionSceneTool::
            createHintFragmentItem(monomer_one_info, monomer_two_info);
        return;
    }
    auto start_qcoords = to_scene_xy(monomer_one_info.pos);
    auto end_qcoords = to_scene_xy(monomer_two_info.pos);
    auto* line_item = new QGraphicsLineItem({start_qcoords, end_qcoords});
    m_invalid_drag_item = line_item;
    line_item->setPen(m_invalid_drag_pen);
    m_scene->addItem(line_item);
    line_item->setVisible(true);
}

void DrawMonomericConnectionSceneTool::addDragStructureToMolModel(
    const HintFragmentMonomerInfo& hint_start_monomer_info,
    const HintFragmentMonomerInfo& hint_end_monomer_info)
{
    if (hint_start_monomer_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG &&
        hint_end_monomer_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG) {
        AbstractDrawMonomerOrMonomericConnectionSceneTool::
            addDragStructureToMolModel(hint_start_monomer_info,
                                       hint_end_monomer_info);
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
                                               NA_BASE_AP_PAIR);
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
