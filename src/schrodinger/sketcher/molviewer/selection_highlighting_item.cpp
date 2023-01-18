#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"

#include <QBrush>
#include <QPen>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

SelectionHighlightingItem::SelectionHighlightingItem() :
    AbstractHighlightingItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setBrush(SELECTION_BACKGROUND_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));
}

} // namespace sketcher
} // namespace schrodinger
