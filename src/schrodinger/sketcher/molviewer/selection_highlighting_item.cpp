#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"

#include <QBrush>
#include <QPen>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

SelectionHighlightingItem::SelectionHighlightingItem() : QGraphicsPathItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setBrush(SELECTION_BACKGROUND_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));

    // nothing is selected yet, so there's no need to paint this
    setVisible(false);
}

void SelectionHighlightingItem::setSelectionPath(const QPainterPath& path)
{
    setPath(path);
    setVisible(!path.isEmpty());
}

} // namespace sketcher
} // namespace schrodinger
