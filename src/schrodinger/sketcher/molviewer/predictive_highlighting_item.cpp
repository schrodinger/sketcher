#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"

#include <QBrush>
#include <QPen>

#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

namespace schrodinger
{
namespace sketcher
{

PredictiveHighlightingItem::PredictiveHighlightingItem() :
    AbstractHighlightingItem()
{
    setPen(PREDICTIVE_HIGHLIGHTING_COLOR);
    setBrush(PREDICTIVE_HIGHLIGHTING_COLOR);
    setZValue(static_cast<qreal>(ZOrder::PREDICTIVE_HIGHLIGHTING));
}

QPainterPath PredictiveHighlightingItem::getPathForItem(
    AbstractGraphicsItem* const item) const
{
    return item->predictiveHighlightingPath();
}

} // namespace sketcher
} // namespace schrodinger
