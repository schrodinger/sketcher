#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include <QPainter>

namespace schrodinger
{
namespace sketcher
{
RotationItem::RotationItem() : QGraphicsPathItem()
{
    setZValue(static_cast<qreal>(ZOrder::ROTATION_HANDLE));
    QPen pen(ROTATION_ITEM_COLOR);
    pen.setWidthF(ROTATION_ITEM_PEN_WIDTH);
    setPen(pen);
    updatePath();
}

void RotationItem::updatePath()
{
    QPainterPath path;
    auto line = getArmLine();
    auto unit_vector = line.unitVector();
    auto dir_vector = unit_vector.p2() - unit_vector.p1();

    path.addEllipse(line.p1(), ROTATION_ITEM_PIVOT_RADIUS,
                    ROTATION_ITEM_PIVOT_RADIUS);
    if (m_draw_rotation_handle) {
        path.moveTo(line.p1() + dir_vector * ROTATION_ITEM_PIVOT_RADIUS);
        path.lineTo(line.p2() - dir_vector * ROTATION_ITEM_HANDLE_RADIUS);
        path.addEllipse(line.p2(), ROTATION_ITEM_HANDLE_RADIUS,
                        ROTATION_ITEM_HANDLE_RADIUS);
    }
    setPath(path);
}

void RotationItem::setPivotPoint(const QPointF& point)
{
    m_pivot_point = point;
    updatePath();
}

QPointF RotationItem::getPivotPoint() const
{
    return m_pivot_point;
}

QLineF RotationItem::getArmLine() const
{
    QLineF line(m_pivot_point, m_pivot_point + QPointF(m_arm_length, 0));
    line.setAngle(m_rotation_angle);
    return line;
}

QPointF RotationItem::getHandlePoint() const
{
    return getArmLine().p2();
}

bool RotationItem::isInsideHandle(const QPointF& point) const
{
    return m_draw_rotation_handle &&
           (QLineF(point, getHandlePoint()).length() <=
            ROTATION_ITEM_HANDLE_RADIUS);
}

void RotationItem::setArmAngle(const float& rotation_angle)
{
    m_rotation_angle = rotation_angle;
    updatePath();
}

float RotationItem::getArmAngle() const
{
    return m_rotation_angle;
}

void RotationItem::setDrawRotationHandle(bool has_handle)
{
    if (has_handle != m_draw_rotation_handle) {
        m_draw_rotation_handle = has_handle;
    }
    updatePath();
}
} // namespace sketcher
} // namespace schrodinger
