#include "schrodinger/sketcher/tool/draw_monomeric_hbond_scene_tool.h"

#include <utility>
#include <variant>

#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

DrawMonomericHBondSceneTool::DrawMonomericHBondSceneTool(const Fonts& fonts,
                                                         Scene* scene,
                                                         MolModel* mol_model) :
    AbstractDrawMonomericConnectionSceneTool(fonts, scene, mol_model)
{
    m_highlight_types = InteractiveItemFlag::MONOMER;
}

QPixmap DrawMonomericHBondSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/monomeric_hbond.svg");
}

void DrawMonomericHBondSceneTool::labelUnboundAttachmentPointsOnMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item,
    QGraphicsItemGroup& attachment_point_labels_group,
    std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
}

std::optional<HintFragmentMonomerInfo>
DrawMonomericHBondSceneTool::getHintFragmentMonomerInfoForDragEnd(
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
        auto [hovered_monomer_item, _drag_end_ap_item] =
            std::get<MonomerAndAPItems>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            hovered_monomer_item, NA_BASE_AP_PAIR);
    }
}

void DrawMonomericHBondSceneTool::updateHoveredUnboundAP(
    UnboundMonomericAttachmentPointItem* hovered_ap_item)
{
}

std::pair<DragEndInfo, AbstractMonomerItem*>
DrawMonomericHBondSceneTool::getDragEndInfo(const QPointF& scene_pos)
{
    auto* hovered_monomer_item = getTopMonomerItemAt(scene_pos);
    if (m_drag_start_monomer_item == nullptr ||
        hovered_monomer_item == nullptr ||
        hovered_monomer_item == m_drag_start_monomer_item ||
        !dragCanFormConnectionTo(hovered_monomer_item, NA_BASE_AP_PAIR)) {
        // we can't form a bond to this location
        return {scene_pos, nullptr};
    } else {
        return {std::make_pair(hovered_monomer_item, nullptr),
                hovered_monomer_item};
    }
}

std::string DrawMonomericHBondSceneTool::getAPModelNameAtPosOverMonomer(
    const QPointF& scene_pos) const
{
    return NA_BASE_AP_PAIR;
}

AbstractMonomerItem*
DrawMonomericHBondSceneTool::getMonomerNear(const QPointF& scene_pos) const
{
    return nullptr;
}

} // namespace sketcher
} // namespace schrodinger
