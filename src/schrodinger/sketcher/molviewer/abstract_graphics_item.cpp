#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

#include <QRectF>

namespace schrodinger
{
namespace sketcher
{

AbstractGraphicsItem::AbstractGraphicsItem(QGraphicsItem* parent) :
    QGraphicsItem(parent)
{
    setFlag(QGraphicsItem::ItemIsSelectable, true);
}

QRectF AbstractGraphicsItem::boundingRect() const
{
    return QRectF(m_bounding_rect);
}

QPainterPath AbstractGraphicsItem::shape() const
{
    return QPainterPath(m_shape);
}

QPainterPath AbstractGraphicsItem::selectionHighlightingPath() const
{
    return QPainterPath(m_selection_highlighting_path);
}

QPainterPath AbstractGraphicsItem::predictiveHighlightingPath() const
{
    return QPainterPath(m_predictive_highlighting_path);
}

} // namespace sketcher
} // namespace schrodinger
