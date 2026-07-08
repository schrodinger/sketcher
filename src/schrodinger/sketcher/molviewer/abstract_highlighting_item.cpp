#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

#include <functional>
#include <vector>

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

/**
 * Calculate the union of all specified paths. This function uses
 * divide-and-conquer to keep the intermediate paths as small as possible, which
 * dramatically improves performance over a naive for loop implementation.
 */
QPainterPath merge_paths(const std::vector<QPainterPath>& paths,
                         const size_t begin, const size_t end)
{
    if (begin == end) {
        return {};
    }
    if (begin + 1 == end) {
        return paths[begin];
    }

    size_t mid = begin + (end - begin) / 2;
    auto left = merge_paths(paths, begin, mid);
    auto right = merge_paths(paths, mid, end);
    return left | right;
}

QPainterPath AbstractHighlightingItem::buildHighlightingPathForItems(
    const QList<const QGraphicsItem*>& items) const
{
    std::vector<QPainterPath> paths;

    for (auto item : items) {
        if (auto* molviewer_item =
                dynamic_cast<const AbstractGraphicsItem*>(item)) {
            auto local_path = getPathForItem(molviewer_item);
            auto scene_path = molviewer_item->mapToScene(local_path);
            paths.push_back(scene_path);
        }
    }
    return merge_paths(paths, 0, paths.size());
}

} // namespace sketcher
} // namespace schrodinger
