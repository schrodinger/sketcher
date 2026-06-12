#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

#include <functional>

#include <QList>

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

void AbstractHighlightingItem::highlightItem(const QGraphicsItem* const item)
{
    highlightItems({item});
}

void AbstractHighlightingItem::highlightItems(
    const QList<QGraphicsItem*>& items)
{
    // Convert to const pointers and delegate to the const version
    QList<const QGraphicsItem*> const_items;
    std::ranges::transform(items, std::back_inserter(const_items),
                           std::identity{});
    highlightItems(const_items);
}

void AbstractHighlightingItem::highlightItems(
    const QList<const QGraphicsItem*>& items)
{
    auto updated_items = updateItemsToHighlight(items);
    QPainterPath path = buildHighlightingPathForItems(updated_items);
    setHighlightingPath(path);
}

QList<const QGraphicsItem*> AbstractHighlightingItem::updateItemsToHighlight(
    const QList<const QGraphicsItem*>& items) const
{
    return items;
}

QPainterPath AbstractHighlightingItem::buildHighlightingPathForItems(
    const QList<const QGraphicsItem*>& items) const
{
    QPainterPath path;
    // Set fill rule to ensure overlapping areas don't create "holes"
    path.setFillRule(Qt::WindingFill);

    for (auto item : items) {
        if (auto* molviewer_item =
                dynamic_cast<const AbstractGraphicsItem*>(item)) {
            QPainterPath local_path = getPathForItem(molviewer_item);
            path.addPath(molviewer_item->mapToScene(local_path));
        }
    }
    return path.simplified();
}

} // namespace sketcher
} // namespace schrodinger
