#pragma once

#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing the selection highlighting in a molviewer
 * Scene.
 */
class SelectionHighlightingItem : public AbstractHighlightingItem
{
  public:
    SelectionHighlightingItem();

    QPainterPath
    getPathForItem(AbstractGraphicsItem* const item) const override;
};

} // namespace sketcher
} // namespace schrodinger
