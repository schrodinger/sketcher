#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

#include <cmath>
#include <memory>

#include <QtMath>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
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
}

std::vector<QGraphicsItem*> DrawMonomerSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_attachment_point_labels_group);
    return items;
}

void DrawMonomerSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MONOMERIC);
    if (item == m_hovered_item) {
        // we've already labeled the item under the cursor, so there's nothing
        // to do
        return;
    }
    clearAttachmentPointsLabels();
    m_hovered_item = item;
    if (item == nullptr) {
        // nothing new to label
        return;
    }

    if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
        // hovering over a monomer
        auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(item);
        const auto* monomer = monomer_item->getAtom();
        labelAttachmentPointsOnMonomer(monomer);
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
        // TODO: mutate monomer
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
    const RDKit::Atom* const monomer)
{
    auto ap_names_and_atoms =
        get_bound_attachment_point_names_and_atoms(monomer);
    for (auto& [ap_name, bound_monomer, is_secondary_connection] :
         ap_names_and_atoms) {
        labelBoundAttachmentPoint(monomer, bound_monomer,
                                  is_secondary_connection, ap_name);
    }
    // TODO: add connector nubs for available attachment points
}

void DrawMonomerSceneTool::labelAttachmentPointsOnConnector(
    const RDKit::Bond* const connector, const bool is_secondary_connection)
{
    auto begin_monomer = connector->getBeginAtom();
    auto end_monomer = connector->getEndAtom();
    auto begin_ap_name = get_attachment_point_name_for_atom(
        begin_monomer, connector, is_secondary_connection);
    auto end_ap_name = get_attachment_point_name_for_atom(
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

/**
 * Position the given rectangle to label a monomer's bound attachment point
 * @param ap_label_rect The rectangle to position. It should already be sized
 * correctly for the attachment point label.
 * @param monomer_coords The coordinates of the monomer being labeled
 * @param bound_coords The coordinates of the other monomer involved in the bond
 */
static void position_ap_label_rect(QRectF& ap_label_rect,
                                   const RDGeom::Point3D& monomer_coords,
                                   const RDGeom::Point3D& bound_coords)
{
    auto qline = QLineF(to_scene_xy(monomer_coords), to_scene_xy(bound_coords));
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

    auto ap_qname = QString::fromStdString(ap_name);
    // convert apostrophes in nucleic acid attachment point names to Unicode
    // primes
    ap_qname.replace('\'', "′");
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
    const int HORIZONTAL_THRESHOLD  = 25;
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
}

} // namespace sketcher
} // namespace schrodinger
