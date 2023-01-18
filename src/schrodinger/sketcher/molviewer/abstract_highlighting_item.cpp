#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

namespace schrodinger
{
namespace sketcher
{

AbstractHighlightingItem::AbstractHighlightingItem() : QGraphicsPathItem()
{
    // nothing is highlighted yet, so there's no need to paint this
    setVisible(false);
}

void AbstractHighlightingItem::setHighlightingPath(const QPainterPath& path)
{
    // note that setPath is a no-op if path == the current path
    setPath(path);
    setVisible(!path.isEmpty());
}

void AbstractHighlightingItem::clearHighlightingPath()
{
    setHighlightingPath(QPainterPath());
}

} // namespace sketcher
} // namespace schrodinger
