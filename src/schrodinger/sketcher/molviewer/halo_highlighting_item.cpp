#include "schrodinger/sketcher/molviewer/halo_highlighting_item.h"

#include <QBrush>
#include <QPen>
#include <QPainter>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

HaloHighlightingItem::HaloHighlightingItem() : AbstractHighlightingItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setBrush(SELECTION_BACKGROUND_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));
}

QPainterPath
HaloHighlightingItem::getPathForItem(AbstractGraphicsItem* const item) const
{
    return item->selectionHighlightingPath();
}

} // namespace sketcher
} // namespace schrodinger
