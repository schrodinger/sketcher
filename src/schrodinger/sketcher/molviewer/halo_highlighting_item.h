#pragma once

#include "schrodinger/sketcher/molviewer/abstract_highlighting_item.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing the colored halo highlighting in a molviewer
 * Scene.
 */
class HaloHighlightingItem : public AbstractHighlightingItem
{
  public:
    HaloHighlightingItem(qreal z_value);

    QPainterPath
    getPathForItem(AbstractGraphicsItem* const item) const override;
};

} // namespace sketcher
} // namespace schrodinger
