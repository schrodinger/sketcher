#pragma once

#include <QGraphicsItem>
#include <QPainterPath>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An abstract parent class for all graphics items that represent the molecule
 * itself (i.e. AtomItem and BondItem, but not SelectionHighlightingItem or
 * PredictiveHighlightingItem) in a molviewer scene.
 */
class SKETCHER_API AbstractGraphicsItem : public QGraphicsItem
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
    QPainterPath shape() const override;

    /**
     * Return the path to draw for selection highlighting when this item is
     * selected.
     */
    QPainterPath selectionHighlightingPath() const;

    /**
     * Return the path to draw for predictive highlighting when this item is
     * hovered over.
     */
    QPainterPath predictiveHighlightingPath() const;

  protected:
    // Type integers for all AbstractGraphicsItem subclasses
    enum class ItemType {
        ATOM = GraphicsItemType::ABSTRACT_GRAPHICS_ITEM_BASE,
        BOND,
        S_GROUP,
        NON_MOLECULAR,
        AMINO_ACID,
        NUCLEIC_ACID_PHOSPHATE,
        NUCLEIC_ACID_SUGAR,
        NUCLEIC_ACID_BASE,
        CHEM_MONOMER,
        MONOMER_CONNECTOR,
    };

    /// The rect to be returned from boundingRect().  Subclasses are responsible
    /// for keeping this value up to date.
    QRectF m_bounding_rect;

    /// The path to be returned from shape().  Subclasses are responsible for
    /// keeping this value up to date.
    QPainterPath m_shape;

    /// The path to be returned from selectionHighlightingPath().  Subclasses
    /// are responsible for keeping this value up to date.
    QPainterPath m_selection_highlighting_path;

    /// The path to be returned from predictiveHighlightingPath().  Subclasses
    /// are responsible for keeping this value up to date.
    QPainterPath m_predictive_highlighting_path;
};

} // namespace sketcher
} // namespace schrodinger
