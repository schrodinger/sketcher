#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

#include <QPainter>

#include "schrodinger/rdkit_extensions/monomer_directions.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Calculate and return the hover area of an attachment point. The hover area is
 * the area that the cursor can click/hover on to draw/hint a connection to the
 * attachment point. It includes the line, a small amount of space on either
 * side of the line, and the label, but excludes any coordinates within the
 * bounding rect of the parent monomer.
 */
static QPainterPath calculate_hover_area(QRectF ap_bounding_rect,
                                         QRectF monomer_bounding_rect,
                                         const Direction ap_direction)
{
    static constexpr int BIG_NUMBER = 500;
    if (ap_direction == Direction::N || ap_direction == Direction::S) {
        // make sure we have at least a little space to the left and right of
        // the attachment point's line
        if (ap_bounding_rect.left() > -UNBOUND_AP_MIN_HOVER_HALF_WIDTH) {
            ap_bounding_rect.setLeft(-UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
        }
        if (ap_bounding_rect.right() < UNBOUND_AP_MIN_HOVER_HALF_WIDTH) {
            ap_bounding_rect.setRight(UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
        }
        // make sure that the hover area doesn't extend above or below the
        // parent monomer. To do this, we extend the monomer's bounding rect in
        // the Y direction. We don't care how far we extend it, so long as
        // monomer_bounding_rect winds up taller than ap_bounding_rect.
        monomer_bounding_rect.adjust(-BIG_NUMBER, 0, BIG_NUMBER, 0);
    } else if (ap_direction == Direction::E || ap_direction == Direction::W) {
        // make sure we have at least a little space above and below the
        // attachment point's line
        if (ap_bounding_rect.top() > -UNBOUND_AP_MIN_HOVER_HALF_WIDTH) {
            ap_bounding_rect.setTop(-UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
        }
        if (ap_bounding_rect.bottom() < UNBOUND_AP_MIN_HOVER_HALF_WIDTH) {
            ap_bounding_rect.setBottom(UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
        }
        // make sure that the hover area doesn't extend to the left or right of
        // the parent monomer. To do this, we extend the monomer's bounding rect
        // in the X direction. We don't care how far we extend it, so long as
        // monomer_bounding_rect winds up wider than ap_bounding_rect.
        monomer_bounding_rect.adjust(0, -BIG_NUMBER, 0, BIG_NUMBER);
    }

    QPainterPath hover_area;
    hover_area.addRect(ap_bounding_rect);
    QPainterPath monomer_bounding_path;
    monomer_bounding_path.addRect(monomer_bounding_rect);
    hover_area -= monomer_bounding_path;
    return hover_area;
}

/**
 * Convert a Direction enum value to a vector in scene coordinates.
 * Note: Y-axis is inverted in Qt (positive Y goes down), so N is (0, -1).
 * @param dir The direction
 * @return vector pointing in that direction
 */
static QPointF direction_to_qt_vector(Direction dir)
{
    auto mol_vec = rdkit_extensions::direction_to_vector(dir);
    return QPointF(mol_vec.x, -mol_vec.y);
}

/**
 * Calculate the geometry required to draw the attachment point
 * @return A tuple of
 *   - coordinates for the end of the attachment point's line
 *   - the text to display for the label
 *   - the rectangle to draw the label in
 *   - the bounding rect of the attachment point
 *   - the hover area for the attachment point  (See the
 *     get_hover_area_for_unbound_monomer_attachment_point_item docstring for an
 *     explanation of hover area)
 */
static std::tuple<QPointF, QString, QRectF, QRectF, QPainterPath>
calculate_geometry(const UnboundAttachmentPoint& attachment_point,
                   const AbstractMonomerItem* const parent_monomer,
                   const Fonts& fonts)
{
    QPointF dir = direction_to_qt_vector(attachment_point.direction);
    QRectF parent_bounds = parent_monomer->boundingRect();
    qreal half_width = parent_bounds.width() / 2.0;
    qreal half_height = parent_bounds.height() / 2.0;

    // Project the direction onto the bounding rect to determine how much of the
    // line will be drawn behind the parent monomer's shape
    qreal extent;
    if (qFuzzyIsNull(dir.x())) {
        // Vertical direction (N or S)
        extent = half_height;
    } else if (qFuzzyIsNull(dir.y())) {
        // Horizontal direction (E or W)
        extent = half_width;
    } else {
        // Diagonal direction - find where the ray exits the bounding rect
        qreal t_x = half_width / qAbs(dir.x());
        qreal t_y = half_height / qAbs(dir.y());
        extent = qMin(t_x, t_y);
    }
    auto line_end = dir * (extent + UNBOUND_AP_LINE_LENGTH);

    auto label_text = prep_attachment_point_name(attachment_point.display_name);
    auto label_rect =
        fonts.m_monomeric_attachment_point_label_fm.boundingRect(label_text);
    position_ap_label_rect(label_rect, {0.0, 0.0}, dir);

    qreal half_line_width = UNBOUND_AP_LINE_THICKNESS / 2.0;
    QRectF line_bounds = QRectF(QPointF(0, 0), line_end).normalized();
    line_bounds.adjust(-half_line_width, -half_line_width, half_line_width,
                       half_line_width);

    qreal radius = UNBOUND_AP_CIRCLE_DIAMETER / 2.0;
    QRectF circle_bounds(line_end.x() - radius, line_end.y() - radius,
                         UNBOUND_AP_CIRCLE_DIAMETER,
                         UNBOUND_AP_CIRCLE_DIAMETER);
    auto bounding_rect = line_bounds.united(circle_bounds).united(label_rect);
    auto hover_area = calculate_hover_area(bounding_rect, parent_bounds,
                                           attachment_point.direction);
    return {line_end, label_text, label_rect, bounding_rect, hover_area};
}

QPainterPath get_hover_area_for_unbound_monomer_attachment_point_item(
    const UnboundAttachmentPoint& attachment_point,
    const AbstractMonomerItem* const parent_monomer, const Fonts& fonts)
{
    auto [line_end, label_text, label_rect, bounding_rect, hover_area] =
        calculate_geometry(attachment_point, parent_monomer, fonts);
    return hover_area;
}

UnboundMonomericAttachmentPointItem::UnboundMonomericAttachmentPointItem(
    const UnboundAttachmentPoint& attachment_point,
    AbstractMonomerItem* parent_monomer, const QColor& color,
    const Fonts& fonts) :
    QGraphicsItem(parent_monomer),
    m_attachment_point(attachment_point),
    m_fonts(&fonts)
{
    setFlag(QGraphicsItem::ItemStacksBehindParent);

    m_line_pen.setWidthF(UNBOUND_AP_LINE_THICKNESS);
    m_line_pen.setCapStyle(Qt::RoundCap);
    m_line_pen.setColor(color);
    m_circle_brush.setColor(color);

    std::tie(m_line_end, m_label_text, m_label_rect, m_bounding_rect,
             m_hover_area) =
        calculate_geometry(m_attachment_point, parent_monomer, *m_fonts);
}

int UnboundMonomericAttachmentPointItem::type() const
{
    return Type;
}

QRectF UnboundMonomericAttachmentPointItem::boundingRect() const
{
    return m_bounding_rect;
}

void UnboundMonomericAttachmentPointItem::paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
    painter->save();

    // Draw the line from center to endpoint
    painter->setPen(m_line_pen);
    painter->drawLine(QPointF(0, 0), m_line_end);

    // Draw the filled circle at the endpoint
    painter->setPen(Qt::NoPen);
    painter->setBrush(m_circle_brush);
    qreal radius = UNBOUND_AP_CIRCLE_DIAMETER / 2.0;
    painter->drawEllipse(m_line_end, radius, radius);

    // Draw the label
    painter->setFont(m_fonts->m_monomeric_attachment_point_label_font);
    painter->setPen(m_line_pen);
    painter->drawText(m_label_rect, Qt::AlignCenter, m_label_text);

    painter->restore();
}

bool UnboundMonomericAttachmentPointItem::withinHoverArea(
    const QPointF& scene_pos) const
{
    return m_hover_area.contains(mapFromScene(scene_pos));
}

const UnboundAttachmentPoint&
UnboundMonomericAttachmentPointItem::getAttachmentPoint() const
{
    return m_attachment_point;
}

} // namespace sketcher
} // namespace schrodinger
