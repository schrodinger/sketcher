#pragma once

#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing the predictive highlighting (i.e. highlighting
 * atoms and bonds as they're moused over) in a molviewer Scene.
 */
class PredictiveHighlightingItem : public AbstractHighlightingItem
{
  public:
    PredictiveHighlightingItem();

    // overriden AbstractHighlightingItem methods
    QPainterPath
    getPathForItem(const AbstractGraphicsItem* const item) const override;

    QList<const QGraphicsItem*> updateItemsToHighlight(
        const QList<const QGraphicsItem*>& items) const override;
};

} // namespace sketcher
} // namespace schrodinger
