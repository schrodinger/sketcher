/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <QGraphicsPathItem>
#include <QPainterPath>

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for drawing the selection highlighting in a molviewer
 * Scene.
 */
class SelectionHighlightingItem : public QGraphicsPathItem
{
  public:
    SelectionHighlightingItem();

    void setSelectionPath(const QPainterPath& path);
};

} // namespace sketcher
} // namespace schrodinger
