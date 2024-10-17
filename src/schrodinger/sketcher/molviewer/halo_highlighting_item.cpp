#include "schrodinger/sketcher/molviewer/halo_highlighting_item.h"

#include <QBrush>
#include <QPen>
#include <QPainter>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

HaloHighlightingItem::HaloHighlightingItem(qreal z_value) :
    AbstractHighlightingItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setBrush(SELECTION_BACKGROUND_COLOR);
    setZValue(z_value);
}

QPainterPath
HaloHighlightingItem::getPathForItem(AbstractGraphicsItem* const item) const
{
    return item->selectionHighlightingPath();
}

} // namespace sketcher
} // namespace schrodinger
