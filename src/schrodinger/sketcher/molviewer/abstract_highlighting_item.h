#pragma once

#include <QGraphicsPathItem>
#include <QPainterPath>

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing highlighting (either selection or predictive)
 * in a molviewer Scene.
 */
class AbstractHighlightingItem : public QGraphicsPathItem
{
  public:
    AbstractHighlightingItem();

    /**
     * Update the path to be highlighted.  If the provided path is empty, then
     * this graphics item will be hidden.  If the path is non-empty, then this
     * graphics item will be shown.
     */
    void setHighlightingPath(const QPainterPath& path);

    /**
     * Clear the path to be highlighted and hide this graphics item.
     */
    void clearHighlightingPath();
};

} // namespace sketcher
} // namespace schrodinger
