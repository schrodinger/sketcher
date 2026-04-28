#include "schrodinger/sketcher/tool/abstract_monomer_scene_tool.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <QGraphicsItem>
#include <QPainterPath>

#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_hint_fragment_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDGeom::Point3D monomer_pos,
                                     const Direction dir)
{
    auto offset = rdkit_extensions::direction_to_vector(dir);
    offset *= BOND_LENGTH;
    return monomer_pos + offset;
}

template <typename T>
void release_after_parent_destroyed(QGraphicsItemGroup& group,
                                    std::vector<T*>& unbound_ap_items)
{
    for (auto* item : group.childItems()) {
        group.removeFromGroup(item);
        delete item;
    }
    unbound_ap_items.clear();
}

// explicitly instanatiate the UnboundMonomericAttachmentPointItem version of
// this function since it's used in DrawMonomerSceneTool
template void release_after_parent_destroyed(
    QGraphicsItemGroup& group,
    std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

AbstractMonomerSceneTool::AbstractMonomerSceneTool(const Fonts& fonts,
                                                   Scene* scene,
                                                   MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_fonts(fonts)
{
    m_attachment_point_labels_group.setZValue(static_cast<qreal>(ZOrder::HINT));
}

AbstractMonomerSceneTool::~AbstractMonomerSceneTool()
{
    // explicitly erase any attachment point labels. Without this, the bound
    // attachment point labels would be deleted regardless when
    // m_attachment_point_labels_group is destroyed, but the unbound
    // attachment point labels are parented to their monomer, so they would
    // outlive the scene tool without this call.
    clearAttachmentPointsLabelsAndHintFragmentItem();
}

std::vector<QGraphicsItem*> AbstractMonomerSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_attachment_point_labels_group);
    return items;
}

void AbstractMonomerSceneTool::updateColorsAfterBackgroundColorChange(
    bool is_dark_mode)
{
    StandardSceneToolBase::updateColorsAfterBackgroundColorChange(is_dark_mode);
    m_unbound_ap_label_color =
        is_dark_mode ? UNBOUND_AP_LABEL_COLOR_DARK_BG : UNBOUND_AP_LABEL_COLOR;
    m_bound_ap_label_color =
        is_dark_mode ? BOUND_AP_LABEL_COLOR_DARK_BG : BOUND_AP_LABEL_COLOR;
}

void AbstractMonomerSceneTool::onStructureUpdated()
{
    release_after_parent_destroyed(m_attachment_point_labels_group,
                                   m_unbound_ap_items);
    m_hovered_item = nullptr;
    m_hovered_ap_item = nullptr;
    clearHintFragmentItem();
    m_predictive_highlighting_item.clearHighlightingPath();
}

bool AbstractMonomerSceneTool::isItemPartOfHintFragment(
    QGraphicsItem* item) const
{
    return m_hint_fragment_item != nullptr &&
           (item == m_hint_fragment_item ||
            item->group() == m_hint_fragment_item);
}

AbstractMonomerItem*
AbstractMonomerSceneTool::getTopMonomerItemAt(const QPointF& scene_pos) const
{
    // check to see if we're over a monomer or unbound attachment point item
    for (auto* item : m_scene->items(scene_pos)) {
        if (isItemPartOfHintFragment(item)) {
            // ignore the fragment hint that was drawn by this class
            continue;
        }
        if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
            return static_cast<AbstractMonomerItem*>(item);
        }
        auto* ap_item =
            qgraphicsitem_cast<UnboundMonomericAttachmentPointItem*>(item);
        if (ap_item && ap_item->withinHoverArea(scene_pos)) {
            return static_cast<AbstractMonomerItem*>(item->parentItem());
        }
    }

    // if we're not over any of those, check to see if we're near a monomer. If
    // we are, check to see whether we'd be over one of its attachment points
    // once they're drawn
    QPainterPath near_scene_pos;
    // We want a circle radius that will definitely include the relevant monomer
    // if we're over an attachment point hover area, but we don't care if it's
    // too large since we'll do further filtering below.  Instead of trying to
    // take the font size, etc. into account, we just pick a number that's close
    // to correct and then double it.
    auto radius =
        2 * std::max(UNBOUND_AP_LINE_LENGTH, UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
    near_scene_pos.addEllipse(scene_pos, radius, radius);
    for (auto* item : m_scene->items(near_scene_pos)) {
        if (!item_matches_type_flag(item, InteractiveItemFlag::MONOMER) ||
            isItemPartOfHintFragment(item)) {
            continue;
        }
        auto* monomer_item = static_cast<AbstractMonomerItem*>(item);
        auto local_pos = monomer_item->mapFromScene(scene_pos);
        auto* monomer = monomer_item->getAtom();
        auto [bound_aps, unbound_aps] =
            get_attachment_points_for_monomer(monomer);
        for (auto cur_unbound_ap : unbound_aps) {
            auto unbound_ap_hover_area =
                get_hover_area_for_unbound_monomer_attachment_point_item(
                    cur_unbound_ap, monomer_item, m_fonts);
            if (unbound_ap_hover_area.contains(local_pos)) {
                return monomer_item;
            }
        }
    }
    return nullptr;
}

MonomerConnectorItem* AbstractMonomerSceneTool::getTopMonomerConnectorItemAt(
    const QPointF& scene_pos) const
{
    for (auto* item : m_scene->items(scene_pos)) {
        if (!isItemPartOfHintFragment(item) &&
            item_matches_type_flag(item,
                                   InteractiveItemFlag::MONOMER_CONNECTOR)) {
            return static_cast<MonomerConnectorItem*>(item);
        }
    }
    return nullptr;
}

void AbstractMonomerSceneTool::drawAttachmentPointLabelsFor(
    QGraphicsItem* const item)
{
    clearAttachmentPointsLabelsAndHintFragmentItem();
    if (item == nullptr) {
        return;
    }
    if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
        // hovering over a monomer
        auto* monomer_item = static_cast<AbstractMonomerItem*>(item);
        const auto* monomer = monomer_item->getAtom();
        labelAttachmentPointsOnHoveredMonomer(monomer, monomer_item);
    } else if (item_matches_type_flag(item,
                                      InteractiveItemFlag::MONOMER_CONNECTOR)) {
        // hovering over a monomeric connector
        auto* connector_item = qgraphicsitem_cast<MonomerConnectorItem*>(item);
        const auto* connector = connector_item->getBond();
        labelAttachmentPointsOnConnector(
            connector, connector_item->isSecondaryConnection());
    }
}

template <typename T>
static void clear_graphics_item_group_and_list(QGraphicsItemGroup& group,
                                               std::vector<T*>& items_list)
{
    for (auto* item : group.childItems()) {
        group.removeFromGroup(item);
        delete item;
    }
    for (auto* item : items_list) {
        delete item;
    }
    items_list.clear();
}

void AbstractMonomerSceneTool::clearAttachmentPointsLabelsAndHintFragmentItem()
{
    clear_graphics_item_group_and_list(m_attachment_point_labels_group,
                                       m_unbound_ap_items);
    m_hovered_ap_item = nullptr;
    clearHintFragmentItem();
}

void AbstractMonomerSceneTool::labelAttachmentPointsOnHoveredMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item)
{
    labelAttachmentPointsOnMonomer(monomer, monomer_item,
                                   m_attachment_point_labels_group,
                                   m_unbound_ap_items);
}

void AbstractMonomerSceneTool::labelAttachmentPointsOnMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item,
    QGraphicsItemGroup& attachment_point_labels_group,
    std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    auto [bound_aps, unbound_aps] = get_attachment_points_for_monomer(monomer);
    for (auto& cur_ap : bound_aps) {
        auto* item = create_label_for_bound_attachment_point(
            monomer, cur_ap.bound_monomer, cur_ap.is_secondary_connection,
            cur_ap.display_name, m_bound_ap_label_color, m_fonts, m_scene);
        if (item != nullptr) {
            attachment_point_labels_group.addToGroup(item);
        }
    }
    for (auto& cur_ap : unbound_aps) {
        auto* item = new UnboundMonomericAttachmentPointItem(
            cur_ap, monomer_item, m_unbound_ap_label_color, m_fonts);
        unbound_ap_items.push_back(item);
    }
}

void AbstractMonomerSceneTool::labelAttachmentPointsOnConnector(
    const RDKit::Bond* const connector, const bool is_secondary_connection)
{
    auto ap_label_items = create_attachment_point_labels_for_connector(
        connector, is_secondary_connection, m_bound_ap_label_color, m_fonts,
        m_scene);
    for (auto* item : ap_label_items) {
        m_attachment_point_labels_group.addToGroup(item);
    }
}

void AbstractMonomerSceneTool::clearHintFragmentItem()
{
    delete m_hint_fragment_item;
    m_hint_fragment_item = nullptr;
}

void AbstractMonomerSceneTool::setCursorHintShown(bool show)
{
    if (show != m_cursor_hint_shown) {
        emit newCursorHintRequested(show ? getDefaultCursorPixmap()
                                         : QPixmap());
        m_cursor_hint_shown = show;
    }
}

} // namespace sketcher
} // namespace schrodinger
