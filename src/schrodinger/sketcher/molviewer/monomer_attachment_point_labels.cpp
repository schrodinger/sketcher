#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"

#include <QGraphicsSimpleTextItem>

#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

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

static QGraphicsItem* create_attachment_point_label(const QString& label,
                                                    const QRectF& label_rect,
                                                    const Fonts& fonts,
                                                    const QColor& color)
{
    auto* label_item = new QGraphicsSimpleTextItem(label);
    label_item->setFont(fonts.m_monomeric_attachment_point_label_font);
    label_item->setBrush({color});
    label_item->setPos(label_rect.topLeft());
    return label_item;
}

QGraphicsItem* create_label_for_bound_attachment_point(
    const RDKit::Atom* const monomer, const RDKit::Atom* const bound_monomer,
    const bool is_secondary_connection, const std::string& ap_name,
    const QColor& color, const Fonts& fonts, const Scene* const scene)
{
    // nothing to do if there's no label (e.g. phosphate attachment points)
    if (ap_name.empty()) {
        return nullptr;
    }

    auto conf = monomer->getOwningMol().getConformer();
    auto monomer_coords = conf.getAtomPos(monomer->getIdx());
    auto bound_coords = conf.getAtomPos(bound_monomer->getIdx());

    auto ap_qname = prep_attachment_point_name(ap_name);
    auto ap_label_rect =
        fonts.m_monomeric_attachment_point_label_fm.boundingRect(ap_qname);
    position_ap_label_rect(ap_label_rect, monomer_coords, bound_coords);
    // if the bond is drawn with an arrowhead at this attachment point (e.g.
    // disulfide bonds or branching monomers), offset the label to account for
    // the arrowhead
    if (attachment_point_is_drawn_with_arrowhead(monomer, bound_monomer,
                                                 is_secondary_connection) &&
        scene != nullptr) {
        // TODO: update this so we don't need the scene, since we can't use it
        //       for hint structures (since it doesn't have the hint's atoms in
        //       its bookkeeping, which is why the hint passes a nullptr for
        //       scene)
        const auto* monomer_item = scene->getGraphicsItemForAtom(monomer);
        auto arrowhead_offset = get_monomer_arrowhead_offset(
            *monomer_item, to_scene_xy(bound_coords));
        ap_label_rect.translate(0, -arrowhead_offset);
    }
    return create_attachment_point_label(ap_qname, ap_label_rect, fonts, color);
}

QGraphicsItem*
create_label_for_center_of_connector(const RDKit::Atom* const begin_monomer,
                                     const RDKit::Atom* const end_monomer,
                                     const QString& label, const QColor& color,
                                     const Fonts& fonts)
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
        fonts.m_monomeric_attachment_point_label_fm.boundingRect(label);
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
    return create_attachment_point_label(label, label_rect, fonts, color);
}

std::vector<QGraphicsItem*> create_attachment_point_labels_for_connector(
    const RDKit::Bond* const connector, const bool is_secondary_connection,
    const QColor& color, const Fonts& fonts, const Scene* const scene)
{
    auto begin_monomer = connector->getBeginAtom();
    auto end_monomer = connector->getEndAtom();
    auto begin_ap_name = get_attachment_point_name_for_connection(
        begin_monomer, connector, is_secondary_connection);
    auto end_ap_name = get_attachment_point_name_for_connection(
        end_monomer, connector, is_secondary_connection);

    std::vector<QGraphicsItem*> ap_label_items;
    if (begin_ap_name == NA_BASE_AP_PAIR && end_ap_name == NA_BASE_AP_PAIR) {
        // for nucleic acid base pairs, only have a single "pair" label since
        // the bond is typically too short to fit two separate labels
        auto* item = create_label_for_center_of_connector(
            begin_monomer, end_monomer, QString::fromStdString(NA_BASE_AP_PAIR),
            color, fonts);
        ap_label_items.push_back(item);
    } else {
        auto* item1 = create_label_for_bound_attachment_point(
            begin_monomer, end_monomer, is_secondary_connection, begin_ap_name,
            color, fonts, scene);
        auto* item2 = create_label_for_bound_attachment_point(
            end_monomer, begin_monomer, is_secondary_connection, end_ap_name,
            color, fonts, scene);
        for (auto* item : {item1, item2}) {
            if (item != nullptr) {
                ap_label_items.push_back(item);
            }
        }
    }
    return ap_label_items;
}

} // namespace sketcher
} // namespace schrodinger
