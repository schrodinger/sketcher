#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

#include <QRectF>

namespace schrodinger
{
namespace sketcher
{

AbstractGraphicsItem::AbstractGraphicsItem(QGraphicsItem* parent) :
    QGraphicsItem(parent)
{
}

QRectF AbstractGraphicsItem::boundingRect() const
{
    return QRectF(m_bounding_rect);
}

} // namespace sketcher
} // namespace schrodinger
