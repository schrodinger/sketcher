#pragma once

#include <QGraphicsPathItem>
#include <QPainterPath>

#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

template <typename T> class QList;

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing highlighting (either selection or predictive)
 * in a molviewer Scene.  Concrete subclasses must provide an implementation for
 * getPathForItem.
 */
class AbstractHighlightingItem : public QGraphicsPathItem
{
  public:
    AbstractHighlightingItem();

    /**
     * Clear any existing highlighting and highlight only the specified graphics
     * item.
     */
    void highlightItem(QGraphicsItem* const item);

    /**
     * Clear any existing highlighting and highlight all specified graphics
     * items.
     */
    void highlightItems(const QList<QGraphicsItem*>& items);

    /**
     * Clear the path to be highlighted and hide this graphics item.
     */
    void clearHighlightingPath();

  protected:
    /**
     * Update the path to be highlighted.  If the provided path is empty, then
     * this graphics item will be hidden.  If the path is non-empty, then this
     * graphics item will be shown.
     */
    void setHighlightingPath(const QPainterPath& path);

    /**
     * Get the painter path to use for highlighting the specified graphics item.
     */
    virtual QPainterPath
    getPathForItem(AbstractGraphicsItem* const item) const = 0;

    /**
     * Modify the list of items to highlight immediately prior to highlighting.
     * The default implementation of this method simply returns the input items,
     * but subclasses may override this method to change which items get
     * highlighted.
     */
    virtual QList<QGraphicsItem*>
    updateItemsToHighlight(const QList<QGraphicsItem*>& items) const;

    /**
     * Create a painter path highlighting all specified graphics items.
     */
    QPainterPath
    buildHighlightingPathForItems(const QList<QGraphicsItem*>& items) const;
};

} // namespace sketcher
} // namespace schrodinger
