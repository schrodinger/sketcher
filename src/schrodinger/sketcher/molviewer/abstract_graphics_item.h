/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <QGraphicsItem>

namespace schrodinger
{
namespace sketcher
{

/**
 * An abstract parent class for all graphics items (i.e. AtomItem and
 * BondItem) in a molviewer scene.
 */
class AbstractGraphicsItem : public QGraphicsItem
{
  public:
    AbstractGraphicsItem(QGraphicsItem* parent = nullptr);

    /**
     * Update all cached data and notify the scene of the changes.  This method
     * must be called whenever the scene updates fonts or settings relevant to
     * this item.
     */
    virtual void updateCachedData() = 0;

    // Overridden QGraphicsItem method
    QRectF boundingRect() const override;

  protected:
    // Type integers for all QGraphicsItem subclasses.  This enum ensures that
    // all graphics items have a different value.
    enum class ItemType {
        // sketcherGraphicalObject::ObjectType uses UserType + 1, so we use +
        // 1000 here to make it more obvious if someone accidentally tries to
        // mix and match old and new graphics items
        ATOM = QGraphicsItem::UserType + 1000,
        BOND,
    };

    QRectF m_bounding_rect;
};

} // namespace sketcher
} // namespace schrodinger
