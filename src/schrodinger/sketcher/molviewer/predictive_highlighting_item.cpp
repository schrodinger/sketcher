#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"

#include <QBrush>
#include <QPen>

#include "schrodinger/sketcher/molviewer/constants.h"

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

} // namespace sketcher
} // namespace schrodinger
