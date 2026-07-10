#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

#include <optional>
#include <utility>
#include <variant>

#include <QGraphicsSceneMouseEvent>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

using rdkit_extensions::Direction;

DrawMonomerSceneTool::DrawMonomerSceneTool(
    const std::string& res_name, const rdkit_extensions::ChainType chain_type,
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    AbstractDrawMonomerOrMonomericConnectionSceneTool(res_name, chain_type,
                                                      fonts, scene, mol_model)
{
}

bool DrawMonomerSceneTool::clickShouldMutate(
    const RDKit::Atom* monomer, const MonomerType monomer_type) const
{
    return (monomer_type == m_monomer_type &&
            get_monomer_res_name(monomer) != m_res_name);
}

void DrawMonomerSceneTool::updateHoveredUnboundAP(
    UnboundMonomericAttachmentPointItem* hovered_ap_item)
{
    // hide the hovered "nubbin" and draw a hint fragment instead
    for (auto* ap_item : m_unbound_ap_items) {
        ap_item->setVisible(hovered_ap_item == nullptr ||
                            ap_item != hovered_ap_item);
    }
    drawBoundMonomerHintFor(hovered_ap_item);
}

void DrawMonomerSceneTool::drawBoundMonomerHintFor(
    UnboundMonomericAttachmentPointItem* const ap_item)
{
    clearHintFragmentItem();
    if (ap_item == nullptr) {
        return;
    }
    const auto* monomer_item =
        static_cast<const AbstractMonomerItem*>(m_hovered_monomeric_item);
    auto hovered_monomer_info =
        createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            monomer_item, ap_item->getAttachmentPoint().model_name);
    auto direction = ap_item->getAttachmentPoint().direction;
    auto new_monomer_info = createHintFragmentMonomerInfoForHintToDirection(
        hovered_monomer_info, direction);
    createHintFragmentItem(hovered_monomer_info, new_monomer_info);
}

void DrawMonomerSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = getTopMonomerItemAt(scene_pos);

    if (item == nullptr) {
        // the click was on empty space, so create a new monomer here
        auto mol_pos = to_mol_xy(scene_pos);
        m_mol_model->addMonomer(m_res_name, m_chain_type, mol_pos);
    } else {
        auto [monomer, monomer_type] = get_monomer_and_type(item);
        std::optional<UnboundAttachmentPoint> clicked_ap;
        auto ap_item = getUnboundAttachmentPointAt(scene_pos, true);
        if (ap_item != nullptr) {
            clicked_ap = ap_item->getAttachmentPoint();
        }

        if (clicked_ap.has_value()) {
            // the user clicked on an attachment point or this monomer has a
            // default attachment point for this tool, so add a new monomer
            // bound to that attachment point
            auto new_monomer_ap_name = get_attachment_point_for_new_monomer(
                monomer_type, clicked_ap->model_name, m_monomer_type);
            auto new_pos = get_default_coords_for_bound_monomer(
                monomer, clicked_ap->direction);
            // the attachment point labels won't be valid once the new monomer
            // is added, so clear them now (otherwise we risk a crash)
            clearAttachmentPointsLabelsAndHintFragmentItem();
            m_mol_model->addBoundMonomer(m_res_name, m_chain_type, new_pos,
                                         new_monomer_ap_name, monomer,
                                         clicked_ap->model_name);

        } else if (clickShouldMutate(monomer, monomer_type)) {
            // the user clicked directly on the monomer and the clicked
            // monomer's residue name is different than the tool's, so we mutate
            // the clicked monomer
            clearAttachmentPointsLabelsAndHintFragmentItem();
            m_mol_model->mutateMonomers({monomer}, m_res_name, m_monomer_type);
        }
    }
}

std::optional<HintFragmentMonomerInfo>
DrawMonomerSceneTool::getHintFragmentMonomerInfoForDragStart()
{
    if (m_drag_start_monomer_item == nullptr) {
        // the drag started over empty space
        if (m_monomer_type == MonomerType::PEPTIDE ||
            m_monomer_type == MonomerType::NA_BASE) {
            return createHintFragmentMonomerInfoForHintFromEmptySpace(
                m_mouse_press_scene_pos);
        } else {
            // it doesn't make much biological sense to connect this monomer
            // type to itself, which is what would typically happen when
            // dragging from empty space, so don't start the drag
            return std::nullopt;
        }
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
DrawMonomerSceneTool::getHintFragmentMonomerInfoForDragEnd(
    const HintFragmentMonomerInfo& hint_start_monomer_info,
    const DragEndInfo& drag_end_info)
{
    if (std::holds_alternative<Direction>(drag_end_info)) {
        auto dir = std::get<Direction>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToDirection(
            hint_start_monomer_info, dir);
    } else {
        auto [hovered_monomer_item, drag_end_ap_item] =
            std::get<MonomerAndAPItems>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            hovered_monomer_item, drag_end_ap_item);
    }
}

DragEndInfo DrawMonomerSceneTool::getDragEndInfoForNonMonomerPos(
    const QPointF& cur_scene_pos) const
{
    const qreal dx = cur_scene_pos.x() - m_mouse_press_scene_pos.x();
    const qreal dy = cur_scene_pos.y() - m_mouse_press_scene_pos.y();
    const qreal abs_dx = std::fabs(dx);
    const qreal abs_dy = std::fabs(dy);
    const qreal sqrt2_plus_1 = std::sqrt(2.0) + 1.0;

    if (abs_dy * sqrt2_plus_1 <= abs_dx) {
        // Within 22.5 degrees of horizontal
        return dx >= 0 ? Direction::E : Direction::W;
    } else if (abs_dy >= abs_dx * sqrt2_plus_1) {
        // Within 22.5 degrees of vertical (Qt +Y is down = South)
        return dy >= 0 ? Direction::S : Direction::N;
    } else {
        // Diagonal
        if (dx >= 0) {
            return dy >= 0 ? Direction::SE : Direction::NE;
        } else {
            return dy >= 0 ? Direction::SW : Direction::NW;
        }
    }
}

QPixmap DrawMonomerSceneTool::createDefaultCursorPixmap() const
{
    // the specific number used here (the "1") doesn't matter - we just need any
    // number to form a proper chain ID
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, 1, false);

    std::shared_ptr<AbstractMonomerItem> monomer_item;
    monomer_item.reset(get_monomer_graphics_item(monomer.get(), *m_fonts));
    monomer_item->setMonomerColors(Qt::GlobalColor::transparent,
                                   CURSOR_HINT_COLOR, CURSOR_HINT_COLOR);
    // make sure that the cursor hint is at least a little smaller than an
    // actual monomer
    auto min_scene_size =
        CURSOR_HINT_IMAGE_SIZE * MONOMER_CURSOR_HINT_MIN_SCENE_SIZE_SCALE;
    return cursor_hint_from_graphics_item(monomer_item.get(), min_scene_size);
}

} // namespace sketcher
} // namespace schrodinger
