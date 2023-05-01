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

    QPainterPath
    getPathForItem(AbstractGraphicsItem* const item) const override;
};

} // namespace sketcher
} // namespace schrodinger
