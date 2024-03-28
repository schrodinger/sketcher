/* -------------------------------------------------------------------------
*  class RotationItem
--------------------------------------------------------------------------- */
#pragma once
#include <QGraphicsItem>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API RotationItem : public QGraphicsPathItem
{
  public:
    RotationItem();

    void updatePath();

    void setPivotPoint(const QPointF& point);
    QPointF getPivotPoint() const;
    void setHandlePoint(const QPointF& point);
    QPointF getHandlePoint() const;

    /*
     * @return whether the given point is inside the handle circle
     */
    bool isInsideHandle(const QPointF& point) const;

    /**
     * Directly assign the angle between the pivot point and the handle.
     *
     * Like on a unit circle, a rotation angle of 0 corresponds to the handle
     * being to the right of the pivot point.
     *
     * @param rotation_angle The rotation angle, in degrees
     */
    void setArmAngle(const float& rotation_angle);

    float getArmAngle() const;

    /**
     * Set whether or not this object should display a handle that the user can
     * grab to rotate the structure. e.g. if a single atom
     * is selected rotation has no meaning and the item is displayed without a
     * handle and only allows translation
     */
    void setDrawRotationHandle(bool has_handle);

  protected:
    QPointF m_pivot_point;
    float m_arm_length = 130;
    float m_rotation_angle = 0;

    /**
     * @return the line from the pivot point to the handle
     */
    QLineF getArmLine() const;

    bool m_draw_rotation_handle = false;
};

} // namespace sketcher
} // namespace schrodinger