#pragma once

#include <QGraphicsEllipseItem>
#include <QGraphicsPathItem>
#include <QGraphicsRectItem>
#include <QPolygonF>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

// A mixin to apply common formatting to all SelectionItem classes
template <typename BASE> class SelectionGraphicsItem : public BASE
{

  public:
    SelectionGraphicsItem(QGraphicsItem* parent = nullptr) : BASE(parent)
    {
        this->setZValue(static_cast<qreal>(ZOrder::DRAG_SELECT_OUTLINE));
    }
};

/**
 * A graphics item used to display the outline of a rectangular marquee
 * selection.
 */
class RectSelectionItem : public SelectionGraphicsItem<QGraphicsRectItem>
{
  public:
    RectSelectionItem(QGraphicsItem* parent = nullptr);
};

/**
 * A graphics item used to display the outline of an elliptical marquee
 * selection.
 */

class EllipseSelectionItem : public SelectionGraphicsItem<QGraphicsEllipseItem>
{
  public:
    EllipseSelectionItem(QGraphicsItem* parent = nullptr);
};

/**
 * A graphics item used to display the outline of a lasso selection.
 *
 * Note that we intentionally inherit from QGraphicsPathItem here instead of
 * QGraphicsPolygonItem because QGraphicsPolygonItem only draws closed polygons,
 * and we want an open polygon.
 */
class LassoSelectionItem : public SelectionGraphicsItem<QGraphicsPathItem>
{
  public:
    LassoSelectionItem(QGraphicsItem* parent = nullptr);

    /**
     * Add a point to the lasso path
     */
    void addPoint(const QPointF& point);

    /**
     * Clear the lasso path
     */
    void clearPath();

  protected:
    QPolygonF m_lasso_polygon;
};

} // namespace sketcher
} // namespace schrodinger