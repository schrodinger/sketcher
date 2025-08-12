#include "schrodinger/sketcher/molviewer/selection_items.h"

#include <QPainterPath>

namespace schrodinger
{
namespace sketcher
{

RectSelectionItem::RectSelectionItem(QGraphicsItem* parent) :
    SelectionGraphicsItem<QGraphicsRectItem>(parent)
{
}

EllipseSelectionItem::EllipseSelectionItem(QGraphicsItem* parent) :
    SelectionGraphicsItem<QGraphicsEllipseItem>(parent)
{
}

LassoSelectionItem::LassoSelectionItem(QGraphicsItem* parent) :
    SelectionGraphicsItem<QGraphicsPathItem>(parent)
{
}

void LassoSelectionItem::addPoint(const QPointF& point)
{
    m_lasso_polygon.append(point);
    QPainterPath path;
    path.addPolygon(m_lasso_polygon);
    setPath(path);
}

void LassoSelectionItem::clearPath()
{
    m_lasso_polygon.clear();
    QPainterPath path;
    setPath(path);
}

} // namespace sketcher
} // namespace schrodinger
