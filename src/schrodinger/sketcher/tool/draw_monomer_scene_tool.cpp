#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

#include <cmath>
#include <memory>

#include <QtMath>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QPainterPath>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

DrawMonomerSceneTool::DrawMonomerSceneTool(
    const std::string& res_name, const rdkit_extensions::ChainType chain_type,
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_res_name(res_name),
    m_chain_type(chain_type),
    m_fonts(fonts)
{
    m_highlight_types = InteractiveItemFlag::MONOMERIC;
    // make sure that the cursor hint font is more easily readable at small
    // size. (Note that m_fonts is a copy of the Scene's fonts, not a reference,
    // so this change won't affect anything else.)
    m_fonts.m_main_label_font.setBold(true);
    m_fonts.updateFontMetrics();
    m_attachment_point_labels_group.setZValue(static_cast<qreal>(ZOrder::HINT));
    if (chain_type == rdkit_extensions::ChainType::PEPTIDE) {
        m_monomer_type = MonomerType::PEPTIDE;
    } else if (chain_type == rdkit_extensions::ChainType::CHEM) {
        m_monomer_type = MonomerType::CHEM;
    } else {
        m_monomer_type = get_na_monomer_type_from_res_name(res_name);
    }
}

DrawMonomerSceneTool::~DrawMonomerSceneTool()
{
    // explicitly erase any attachment point labels. Without this, the bound
    // attachment point labels would be deleted regardless when
    // m_attachment_point_labels_group is destroyed, but the unbound attachment
    // point labels are parented to their monomer, so they would outlive the
    // scene tool without this call.
    clearAttachmentPointsLabels();
}

std::vector<QGraphicsItem*> DrawMonomerSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_attachment_point_labels_group);
    return items;
}

QGraphicsItem*
DrawMonomerSceneTool::getTopMonomericItemAt(const QPointF& scene_pos)
{
    // check to see if we're over a monomer, monomeric connector, or unbound
    // attachment point item
    for (auto* item : m_scene->items(scene_pos)) {
        if (item_matches_type_flag(item, InteractiveItemFlag::MONOMERIC)) {
            return item;
        } else if (auto* ap_item =
                       qgraphicsitem_cast<UnboundMonomericAttachmentPointItem*>(
                           item)) {
            if (ap_item->withinHoverArea(scene_pos)) {
                return item->parentItem();
            }
        }
    }

    // if we're not over any of those, check to see if we're near a monomer. If
    // we are, check to see whether we'd be over one of its attachment points
    // once they're drawn
    QPainterPath near_scene_pos;
    // the attachment point label can stick out past the attachment point line,
    // so make the circle a bit bigger than just the line length
    near_scene_pos.addEllipse(scene_pos, 2 * UNBOUND_AP_LINE_LENGTH,
                              2 * UNBOUND_AP_LINE_LENGTH);
    for (auto* item : m_scene->items(near_scene_pos)) {
        if (!item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
            continue;
        }
        const auto* monomer_item =
            dynamic_cast<const AbstractMonomerItem*>(item);
        auto local_pos = monomer_item->mapFromScene(scene_pos);
        auto* monomer = monomer_item->getAtom();
        auto [bound_aps, unbound_aps] =
            get_attachment_points_for_monomer(monomer);
        for (auto cur_unbound_ap : unbound_aps) {
            auto unbound_ap_bounding_rect =
                get_bounding_rect_for_unbound_monomer_attachment_point_item(
                    cur_unbound_ap, monomer_item, m_fonts);
            if (unbound_ap_bounding_rect.contains(local_pos)) {
                return item;
            }
        }
    }

    return nullptr;
}

/**
 * @return the unbound attachment point graphics item representing the
 * attachment point with the "best" number. "Best" is defined using the
 * preferred_order list, with earlier numbers in the list preferred over later
 * numbers. Will return nullptr if no attachment points on the preferred_order
 * list are found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_preferred_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::vector<int>& preferred_order)
{
    auto index_of = [&preferred_order](const auto* ap_item) {
        auto it = std::find(preferred_order.begin(), preferred_order.end(),
                            ap_item->getAttachmentPoint().num);
        auto dist = std::distance(preferred_order.begin(), it);
        // we know that dist is positive
        return static_cast<std::size_t>(dist);
    };
    auto min_it = std::min_element(
        unbound_ap_items.begin(), unbound_ap_items.end(),
        [&index_of](const auto* ap_item_left, const auto* ap_item_right) {
            return index_of(ap_item_left) < index_of(ap_item_right);
        });
    // make sure that the attachment point we found is actually on the
    // preferred_order list
    if (index_of(*min_it) < preferred_order.size()) {
        return *min_it;
    }
    return nullptr;
}

/**
 * Return the unbound attachment point graphics item that represents the
 * attachment point with the lowest number. An attachment points with a custom
 * name will be returned only if there are no numbered attachment point.
 * @param unbound_ap_items The non-empty list of unbound attachment point
 * graphics item
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_min_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    auto min_it = std::min_element(
        unbound_ap_items.begin(), unbound_ap_items.end(),
        [](const auto* ap_item_left, const auto* ap_item_right) {
            auto left_num = ap_item_left->getAttachmentPoint().num;
            auto right_num = ap_item_right->getAttachmentPoint().num;
            return right_num < 0 || left_num < right_num;
        });
    return *min_it;
}

/**
 * @return the unbound attachment point graphics item representing an attachment
 * point with the given name. Will return nullptr if no such attachment point is
 * found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_attachment_point_with_name(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::string& name)
{
    auto it =
        std::find_if(unbound_ap_items.begin(), unbound_ap_items.end(),
                     [&name](const auto* ap_item) {
                         return (ap_item->getAttachmentPoint().name == name);
                     });
    if (it == unbound_ap_items.end()) {
        return nullptr;
    }
    return *it;
}

/**
 * @return the unbound attachment point graphics item representing an attachment
 * point with the given number. Will return nullptr if no such attachment point
 * is found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_attachment_point_with_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const int num)
{
    auto it =
        std::find_if(unbound_ap_items.begin(), unbound_ap_items.end(),
                     [&num](const auto* ap_item) {
                         return (ap_item->getAttachmentPoint().num == num);
                     });
    if (it == unbound_ap_items.end()) {
        return nullptr;
    }
    return *it;
}

UnboundMonomericAttachmentPointItem* get_default_attachment_point(
    const MonomerType hovered_type, const MonomerType tool_type,
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    if (unbound_ap_items.empty()) {
        return nullptr;
    } else if (hovered_type == MonomerType::CHEM) {
        return find_min_attachment_point_by_num(unbound_ap_items);
    } else if (hovered_type == MonomerType::PEPTIDE) {
        if (tool_type == MonomerType::PEPTIDE) {
            return find_preferred_attachment_point_by_num(unbound_ap_items,
                                                          {2, 1, 3});
        } else if (tool_type == MonomerType::CHEM) {
            return find_attachment_point_with_num(unbound_ap_items, 3);
        }
    } else if (hovered_type == MonomerType::NA_BASE) {
        if (tool_type == MonomerType::NA_BASE ||
            tool_type == MonomerType::CHEM) {
            return find_attachment_point_with_name(unbound_ap_items, "pair");
        } else if (tool_type == MonomerType::NA_SUGAR) {
            return find_attachment_point_with_num(unbound_ap_items, 1);
        }
    } else if (hovered_type == MonomerType::NA_SUGAR) {
        if (tool_type == MonomerType::NA_BASE) {
            return find_attachment_point_with_num(unbound_ap_items, 3);
        } else if (tool_type == MonomerType::NA_PHOSPHATE) {
            return find_preferred_attachment_point_by_num(unbound_ap_items,
                                                          {2, 1});
        }
    } else if (hovered_type == MonomerType::NA_PHOSPHATE) {
        if (tool_type == MonomerType::NA_SUGAR) {
            return find_preferred_attachment_point_by_num(unbound_ap_items,
                                                          {2, 1});
        }
    }
    return nullptr;
}

UnboundMonomericAttachmentPointItem*
DrawMonomerSceneTool::getUnboundAttachmentPointAt(const QPointF& scene_pos)
{
    for (auto* ap_item : m_unbound_ap_items) {
        if (ap_item->withinHoverArea(scene_pos)) {
            return ap_item;
        }
    }
    const auto* monomer_item =
        static_cast<const AbstractMonomerItem*>(m_hovered_item);
    auto* monomer = monomer_item->getAtom();
    auto monomer_type = get_monomer_type(monomer);
    if (monomer_type == m_monomer_type &&
        get_monomer_res_name(monomer) != m_res_name) {
        // a click on this monomer should mutate the residue type of the
        // monomer, not add a connection, so we don't select an attachment point
        // when the monomer itself is hovered
        return nullptr;
    }
    return get_default_attachment_point(monomer_type, m_monomer_type,
                                        m_unbound_ap_items);
}

void DrawMonomerSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto* item = getTopMonomericItemAt(scene_pos);
    if (item != m_hovered_item) {
        m_hovered_item = item;
        drawAttachmentPointLabelsFor(item);
    }

    if (!m_unbound_ap_items.empty()) {
        // if we're over a monomer with attachment points, update which
        // attachment point is hovered
        auto* active_ap_item = getUnboundAttachmentPointAt(scene_pos);
        for (auto* ap_item : m_unbound_ap_items) {
            ap_item->setActive(ap_item == active_ap_item);
        }
    }
}

void DrawMonomerSceneTool::drawAttachmentPointLabelsFor(
    QGraphicsItem* const item)
{
    clearAttachmentPointsLabels();
    if (item == nullptr) {
        return;
    }
    if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
        // hovering over a monomer
        auto* monomer_item = static_cast<AbstractMonomerItem*>(item);
        const auto* monomer = monomer_item->getAtom();
        labelAttachmentPointsOnMonomer(monomer, monomer_item);
    } else {
        // hovering over a monomeric connector
        auto* connector_item = qgraphicsitem_cast<MonomerConnectorItem*>(item);
        const auto* connector = connector_item->getBond();
        labelAttachmentPointsOnConnector(
            connector, connector_item->isSecondaryConnection());
    }
}

void DrawMonomerSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MONOMERIC);
    if (item == nullptr) {
        // the click was on empty space, so create a new monomer here
        auto mol_pos = to_mol_xy(scene_pos);
        m_mol_model->addMonomer(m_res_name, m_chain_type, mol_pos);
    } else if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
        // TODO: mutate monomer, or add a new monomer with a connection
        // auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(item);
        // const auto* monomer = monomer_item->getAtom();
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
    monomer_item.reset(get_monomer_graphics_item(monomer.get(), m_fonts));
    monomer_item->setMonomerColors(Qt::GlobalColor::transparent,
                                   CURSOR_HINT_COLOR, CURSOR_HINT_COLOR);
    // make sure that the cursor hint is at least a little smaller than an
    // actual monomer
    auto min_scene_size =
        CURSOR_HINT_IMAGE_SIZE * MONOMER_CURSOR_HINT_MIN_SCENE_SIZE_SCALE;
    return cursor_hint_from_graphics_item(monomer_item.get(), min_scene_size);
}

void DrawMonomerSceneTool::labelAttachmentPointsOnMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item)
{
    auto [bound_aps, unbound_aps] = get_attachment_points_for_monomer(monomer);
    for (auto& cur_ap : bound_aps) {
        labelBoundAttachmentPoint(monomer, cur_ap.bound_monomer,
                                  cur_ap.is_secondary_connection, cur_ap.name);
    }
    for (auto& cur_ap : unbound_aps) {
        auto* item = new UnboundMonomericAttachmentPointItem(
            cur_ap, monomer_item, m_fonts);
        m_unbound_ap_items.push_back(item);
    }
}

void DrawMonomerSceneTool::labelAttachmentPointsOnConnector(
    const RDKit::Bond* const connector, const bool is_secondary_connection)
{
    auto begin_monomer = connector->getBeginAtom();
    auto end_monomer = connector->getEndAtom();
    auto begin_ap_name = get_attachment_point_name_for_connection(
        begin_monomer, connector, is_secondary_connection);
    auto end_ap_name = get_attachment_point_name_for_connection(
        end_monomer, connector, is_secondary_connection);

    if (begin_ap_name == "pair" && end_ap_name == "pair") {
        // for nucleic acid base pairs, only have a single "pair" label since
        // the bond is typically too short to fit two separate labels
        labelCenterOfConnector(begin_monomer, end_monomer, "pair");
    } else {
        labelBoundAttachmentPoint(begin_monomer, end_monomer,
                                  is_secondary_connection, begin_ap_name);
        labelBoundAttachmentPoint(end_monomer, begin_monomer,
                                  is_secondary_connection, end_ap_name);
    }
}

void position_ap_label_rect(QRectF& ap_label_rect,
                            const QPointF& monomer_coords,
                            const QPointF& bound_coords)
{
    auto qline = QLineF(monomer_coords, bound_coords);
    auto angle = qline.angle();

    bool rotate_ccw;
    bool left_aligned;
    bool bottom_aligned;

    // Qt measures angles in [0, 360) degrees, counter-clockwise from a point on
    // the x-axis to the right of the origin.
    if (angle == 0 || angle >= 315) {
        rotate_ccw = true;
        left_aligned = true;
        bottom_aligned = true;
    } else if (angle < 45) {
        rotate_ccw = false;
        left_aligned = true;
        bottom_aligned = false;
    } else if (angle < 90) {
        rotate_ccw = true;
        left_aligned = false;
        bottom_aligned = true;
    } else if (angle <= 135) {
        rotate_ccw = false;
        left_aligned = true;
        bottom_aligned = true;
    } else if (angle < 180) {
        rotate_ccw = true;
        left_aligned = false;
        bottom_aligned = false;
    } else if (angle <= 225) {
        rotate_ccw = false;
        left_aligned = false;
        bottom_aligned = true;
    } else if (angle <= 270) {
        rotate_ccw = true;
        left_aligned = true;
        bottom_aligned = false;
    } else {
        rotate_ccw = false;
        left_aligned = false;
        bottom_aligned = false;
    }

    if (rotate_ccw) {
        angle += MONOMERIC_ATTACHMENT_POINT_LABEL_ANGLE_OFFSET;
    } else {
        angle -= MONOMERIC_ATTACHMENT_POINT_LABEL_ANGLE_OFFSET;
    }
    if (left_aligned) {
        ap_label_rect.moveLeft(0.0);
    } else {
        ap_label_rect.moveRight(0.0);
    }
    if (bottom_aligned) {
        ap_label_rect.moveBottom(0.0);
    } else {
        ap_label_rect.moveTop(0.0);
    }

    qline.setAngle(angle);
    // TODO: explicitly account for monomer size
    qline.setLength(MONOMERIC_ATTACHMENT_POINT_LABEL_DIST);
    ap_label_rect.translate(qline.p2());
}

/**
 * @overload Accepts RDKit coordinates instead of Scene coordinates
 */
static void position_ap_label_rect(QRectF& ap_label_rect,
                                   const RDGeom::Point3D& monomer_coords,
                                   const RDGeom::Point3D& bound_coords)
{
    position_ap_label_rect(ap_label_rect, to_scene_xy(monomer_coords),
                           to_scene_xy(bound_coords));
}

/**
 * For the bond between monomer and bound_monomer, return whether monomer's
 * attachment point has been drawn with an arrowhead (which happens for
 * disulfide bonds and branch monomers).
 */
static bool
attachment_point_is_drawn_with_arrowhead(const RDKit::Atom* const monomer,
                                         const RDKit::Atom* const bound_monomer,
                                         const bool is_secondary_connection)
{
    const auto* bond = monomer->getOwningMol().getBondBetweenAtoms(
        monomer->getIdx(), bound_monomer->getIdx());
    auto [begin_arrowhead, end_arrowhead] =
        does_connector_have_arrowheads(bond, is_secondary_connection);
    if (bond->getBeginAtom() == monomer) {
        return begin_arrowhead;
    } else {
        return end_arrowhead;
    }
}

QString prep_attachment_point_name(const std::string& name)
{
    auto qname = QString::fromStdString(name);
    // convert apostrophes in nucleic acid attachment point names to Unicode
    // primes
    qname.replace('\'', "′");
    return qname;
}

void DrawMonomerSceneTool::labelBoundAttachmentPoint(
    const RDKit::Atom* const monomer, const RDKit::Atom* const bound_monomer,
    const bool is_secondary_connection, const std::string& ap_name)
{
    // nothing to do if there's no label (e.g. phosphate attachment points)
    if (ap_name.empty()) {
        return;
    }

    auto conf = monomer->getOwningMol().getConformer();
    auto monomer_coords = conf.getAtomPos(monomer->getIdx());
    auto bound_coords = conf.getAtomPos(bound_monomer->getIdx());

    auto ap_qname = prep_attachment_point_name(ap_name);
    auto ap_label_rect =
        m_fonts.m_monomeric_attachment_point_label_fm.boundingRect(ap_qname);
    position_ap_label_rect(ap_label_rect, monomer_coords, bound_coords);
    // if the bond is drawn with an arrowhead at this attachment point (e.g.
    // disulfide bonds or branching monomers), offset the label to account for
    // the arrowhead
    if (attachment_point_is_drawn_with_arrowhead(monomer, bound_monomer,
                                                 is_secondary_connection)) {
        const auto* monomer_item = m_scene->getGraphicsItemForAtom(monomer);
        auto arrowhead_offset = get_monomer_arrowhead_offset(
            *monomer_item, to_scene_xy(bound_coords));
        ap_label_rect.translate(0, -arrowhead_offset);
    }
    addAttachmentPointLabel(ap_qname, ap_label_rect);
}

void DrawMonomerSceneTool::addAttachmentPointLabel(const QString& label,
                                                   const QRectF& label_rect)
{
    auto* label_item = new QGraphicsSimpleTextItem(label);
    label_item->setFont(m_fonts.m_monomeric_attachment_point_label_font);
    label_item->setBrush({Qt::black});
    label_item->setPos(label_rect.topLeft());
    m_attachment_point_labels_group.addToGroup(label_item);
}

void DrawMonomerSceneTool::labelCenterOfConnector(
    const RDKit::Atom* const begin_monomer,
    const RDKit::Atom* const end_monomer, const QString& label)
{
    auto conf = begin_monomer->getOwningMol().getConformer();
    auto begin_coords = conf.getAtomPos(begin_monomer->getIdx());
    auto begin_scene_coords = to_scene_xy(begin_coords);
    auto end_coords = conf.getAtomPos(end_monomer->getIdx());
    auto end_scene_coords = to_scene_xy(end_coords);
    auto connector_line = QLineF(begin_scene_coords, end_scene_coords);
    auto connector_angle = connector_line.angle();
    // ensure that the angle of the connector is 0-180, as this will simplify
    // the logic for positioning the label
    if (connector_angle >= 180) {
        // reverse the line if the angle is 180-360
        connector_line = QLineF(end_scene_coords, begin_scene_coords);
        connector_angle = connector_line.angle();
    }

    // we'll place the label just off of the connector's midpoint
    auto midpoint = connector_line.pointAt(0.5);
    QLineF perpendicular = connector_line.normalVector();
    perpendicular.translate(midpoint - perpendicular.p1());
    // make sure that the perpendicular extends above the line, not below. For a
    // vertical line, extend it to the right
    perpendicular.setLength((connector_angle < 90 ? 1 : -1) *
                            MONOMERIC_PAIR_CONNECTOR_LABEL_DIST);
    auto label_rect =
        m_fonts.m_monomeric_attachment_point_label_fm.boundingRect(label);
    QPointF label_pos = perpendicular.p2();
    label_rect.moveCenter(label_pos);

    // a connector within HORIZONTAL_THRESHOLD degrees of the X axis will be
    // treated as horizontal
    const int HORIZONTAL_THRESHOLD = 25;
    bool is_horizontal = connector_angle <= HORIZONTAL_THRESHOLD ||
                         connector_angle >= (180 - HORIZONTAL_THRESHOLD);
    if (is_horizontal) {
        // center the bottom edge on label_pos
        label_rect.moveBottom(label_pos.y());
    } else if (connector_angle < 90) {
        // the connector goes from the lower left to the upper right, so center
        // the right edge on label_pos
        label_rect.moveRight(label_pos.x());
    } else {
        // the connector goes from the lower right to the upper left, so center
        // the left edge on label_pos
        label_rect.moveLeft(label_pos.x());
    }
    addAttachmentPointLabel(label, label_rect);
}

void DrawMonomerSceneTool::clearAttachmentPointsLabels()
{
    for (auto* item : m_attachment_point_labels_group.childItems()) {
        m_attachment_point_labels_group.removeFromGroup(item);
        delete item;
    }
    for (auto* item : m_unbound_ap_items) {
        delete item;
    }
    m_unbound_ap_items.clear();
}

} // namespace sketcher
} // namespace schrodinger
