/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <vector>

#include <QGraphicsItem>
#include <QLineF>
#include <QPen>
#include <QRectF>

#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"

class QPainter;
class QPainterPath;
class QStyleOptionGraphicsItem;
class QWidget;

namespace RDKit
{
class Bond;
}

namespace schrodinger
{
namespace sketcher
{

class AtomItem;

/**
 * A collection of lines and polygons that should be painted using the specified
 * pen.  The pen should be either BondItem::m_solid_pen or
 * BondItem::m_dashed_pen.
 */
struct ToPaint {
    // the pen must be const in order for nested initializer lists (i.e. the
    // syntax used in calculateLinesToPaint) to work
    const QPen& pen;
    std::vector<QLineF> lines;
    std::vector<QPolygon> polygons;
};

/**
 * A Qt graphics item for representing bonds in a molviewer Scene.
 */
class BondItem : public AbstractGraphicsItem
{
  public:
    /**
     * Note that this class does not take ownership of bond, settings, or the
     * atom items.  These objects must not be destroyed while the AtomItem is in
     * use.
     *
     * @param bond The RDKit bond that this item should represent.
     *
     * @param start_atom_item The atom item representing one of the bound atoms.
     *
     * @param end_atom_item The atom item representing the other bound atom.
     *
     * @param settings The settings to be used for painting this item.  To
     * change settings, the calling code (i.e. the Scene) must update this
     * object and then call updateCachedData.
     *
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     * @pre bond != nullptr
     * @pre bond->hasOwningMol()
     */
    BondItem(RDKit::Bond* bond, const AtomItem& start_atom_item,
             const AtomItem& end_atom_item, BondItemSettings& settings,
             QGraphicsItem* parent = nullptr);

    // Type and type() are required for qgraphicsitem_cast.  Note that this
    // implementation is based on the sample code from the QGraphicsItem::Type
    // documentation.
    enum { Type = static_cast<int>(ItemType::BOND) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;
    QPainterPath shape() const override;

  protected:
    // Creating a shared_ptr to an RDKit Bond (or Atom) implicitly creates a
    // copy of the Bond, which means that the new Bond is no longer part of the
    // original molecule, which leads to problems.  Because of this, we store a
    // raw pointer instead.  The RDKit molecule takes care of the lifetime of
    // the Bond, so a BondItem instance must be deleted as soon as its
    // associated Bond is deleted.
    //
    // Also note that m_bond, m_start_item, and m_end_item should only be
    // accessed from within updateCachedData to ensure that we can properly
    // notify the scene of any BondItem changes *before* they happen.
    RDKit::Bond* const m_bond;
    const AtomItem& m_start_item;
    const AtomItem& m_end_item;
    std::vector<ToPaint> m_to_paint;
    QPen m_solid_pen;
    QPen m_dashed_pen;
    QPainterPath m_shape;

    BondItemSettings& m_settings;

    /**
     * Calculate the lines and polygons needed to paint this bond.  Note that
     * this method does *not* do any actual painting.  (The output of this
     * method should be stored in m_to_paint, and the actual painting will be
     * done in `paint`.)
     *
     * @param bond_line A line from m_start_item to m_end_item
     * @return All of the lines and polygons to paint, along with the pens to be
     * used for painting.
     */
    std::vector<ToPaint> calculateLinesToPaint(const QLineF& bond_line) const;

    /**
     * Calculate the lines needed to paint a triple bond.
     *
     * @param trimmed_line A line from m_start_item to m_end_item that has been
     * trimmed so it doesn't overlap the atom item labels
     * @return The lines to paint.  This list will always contain exacatly three
     * items.
     */
    std::vector<QLineF> calcTripleBondLines(const QLineF& trimmed_line) const;

    /**
     * Trim a line so that it doesn't overlap the atom item labels from either
     * bound atom.
     *
     * @param line[in,out] The line to trim
     */
    void trimLineToBoundAtoms(QLineF& line) const;

    /**
     * Trim a line so that it ends at least 4 pixels outside of the given
     * rectangle.  Note that this method assumes that the line is either
     * completely outside of the rectangle (in which case no trimming is
     * required), or that the line has exactly one endpoint within the rectangle
     * (in which case that endpoint is moved outside of the rectangle).  We
     * assume that the line is *not* completely contained within the rectangle
     * and that the line does not pass through the rectangle and come out the
     * other side.
     *
     * @param line[in,out] The line to trim
     * @param subrect[in] The rectangle to trim to
     */
    void trimLineToRect(QLineF& line, const QRectF& subrect) const;

    /**
     * Calculate the intersection point (if any) of a line and a rectangle.  See
     * caveats in the `trimLineToRect` docstring.
     *
     * @param line[in] The line
     * @param rect[in] The rectangle
     * @param i[out] Will be set to the intersection point, if there is one
     * @return Whether the line intersects the rectangle
     */
    bool intersectionOfLineAndRect(const QLineF& line, const QRectF& rect,
                                   QPointF& i) const;

    /**
     * Calculate a painter path around the perimeter of `line` if `line` were
     * drawn with a pen width of `2 * half_width`.
     *
     * @param line The line for the path to go around
     * @param half_width The desired half-width for the painter path
     */
    QPainterPath pathAroundLine(const QLineF& line,
                                const qreal half_width) const;
};

} // namespace sketcher
} // namespace schrodinger
