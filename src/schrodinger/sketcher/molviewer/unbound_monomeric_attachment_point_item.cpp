#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

#include <cmath>

#include <QPainter>

#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Convert a Direction enum value to a unit vector in scene coordinates.
 * Note: Y-axis is inverted in Qt (positive Y goes down), so N is (0, -1).
 * @param dir The direction
 * @return Unit vector pointing in that direction
 */
static QPointF direction_to_unit_vector(Direction dir)
{
    switch (dir) {
        case Direction::N:
            return QPointF(0, -1);
        case Direction::S:
            return QPointF(0, 1);
        case Direction::E:
            return QPointF(1, 0);
        case Direction::W:
            return QPointF(-1, 0);
        case Direction::NE:
            return QPointF(M_SQRT1_2, -M_SQRT1_2);
        case Direction::NW:
            return QPointF(-M_SQRT1_2, -M_SQRT1_2);
        case Direction::SE:
            return QPointF(M_SQRT1_2, M_SQRT1_2);
        case Direction::SW:
            return QPointF(-M_SQRT1_2, M_SQRT1_2);
        default:
            return QPointF(1, 0);
    }
}

static std::tuple<QPointF, QString, QRectF, QRectF>
calculate_geometry(const UnboundAttachmentPoint& attachment_point,
                   const AbstractMonomerItem* const parent_monomer,
                   const Fonts& fonts)
{
    QPointF dir = direction_to_unit_vector(attachment_point.direction);
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

    auto label_text = prep_attachment_point_name(attachment_point.name);
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

    return {line_end, label_text, label_rect, bounding_rect};
}

QRectF get_bounding_rect_for_unbound_monomer_attachment_point_item(
    const UnboundAttachmentPoint& attachment_point,
    const AbstractMonomerItem* const parent_monomer, const Fonts& fonts)
{
    auto [line_end, label_text, label_rect, bounding_rect] =
        calculate_geometry(attachment_point, parent_monomer, fonts);
    return bounding_rect;
}

UnboundMonomericAttachmentPointItem::UnboundMonomericAttachmentPointItem(
    const UnboundAttachmentPoint& attachment_point,
    AbstractMonomerItem* parent_monomer, const Fonts& fonts) :
    QGraphicsItem(parent_monomer),
    m_attachment_point(attachment_point),
    m_fonts(&fonts)
{
    setFlag(QGraphicsItem::ItemStacksBehindParent);

    m_line_pen.setWidthF(UNBOUND_AP_LINE_THICKNESS);
    m_line_pen.setCapStyle(Qt::RoundCap);

    std::tie(m_line_end, m_label_text, m_label_rect, m_bounding_rect) =
        calculate_geometry(m_attachment_point, parent_monomer, *m_fonts);
    m_hover_area.addRect(m_bounding_rect);
    QPainterPath parent_bounds_path;
    parent_bounds_path.addRect(parent_monomer->boundingRect());
    m_hover_area -= mapFromParent(parent_bounds_path);
    updateColors();
}

int UnboundMonomericAttachmentPointItem::type() const
{
    return Type;
}

void UnboundMonomericAttachmentPointItem::setActive(bool active)
{
    if (m_is_active != active) {
        m_is_active = active;
        updateColors();
    }
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

void UnboundMonomericAttachmentPointItem::updateColors()
{
    QColor color =
        m_is_active ? UNBOUND_AP_ACTIVE_COLOR : UNBOUND_AP_INACTIVE_COLOR;
    m_line_pen.setColor(color);
    m_circle_brush.setColor(color);
    update();
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
