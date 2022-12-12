/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <vector>
#include <tuple>
#include <utility>

#include <QtGlobal>
#include <QBrush>
#include <QGraphicsItem>
#include <QLineF>
#include <QPen>
#include <QRectF>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"

class QPainter;
class QPainterPath;
class QStyleOptionGraphicsItem;
class QWidget;

namespace RDKit
{
class Atom;
class Bond;
class Conformer;
class RingInfo;
class ROMol;
} // namespace RDKit

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
    std::vector<QPolygonF> polygons;
};

/**
 * A Qt graphics item for representing bonds in a molviewer Scene.
 */
class SKETCHER_API BondItem : public AbstractGraphicsItem
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
    QBrush m_solid_brush = QBrush(Qt::black);
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
     * Calculate the polygon for an arrow tip on the end of trimmed_line
     *
     * @param trimmed_line The line to draw the arrow tip on
     * @return The newly generated polygon
     */
    QPolygonF calcArrowTip(const QLineF& trimmed_line) const;

    /**
     * Calculate the three endpoints of the wedge for a a chiral bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return A tuple of endpoints.  The first point is the narrow end of the
     * wedge.
     */
    std::tuple<QPointF, QPointF, QPointF>
    calcWedgeEnds(const QLineF& trimmed_line) const;

    /**
     * Calculate the polygon for a solid wedge representing an "up" chiral bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return The requested polygon
     */
    QPolygonF calcSolidWedge(const QLineF& trimmed_line) const;

    /**
     * Calculate all of the values needed to generate a dashed wedge bond (a
     * "down" chiral bond) or a squiggly wedge bond (an "up or down" chiral
     * bond).
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return A tuple of
     *   - the number of dashes to draw in the wedge
     *   - the point at the narrow end of the wedge
     *   - the direction from the narrow end of the wedge to one side of the
     *     wide end of the wedge
     *   - the direction from the narrow end of the wedge to the other side of
     *     the wide end of the wedge
     *   - a vector to store the generated dashes in.  This vector is empty but
     *     memory has been reserved appropriately.
     */
    std::tuple<unsigned int, QPointF, QPointF, QPointF, std::vector<QLineF>>
    calcDashedWedgeParams(const QLineF& trimmed_line) const;

    /**
     * Calculate the lines for a dashed wedge representing a "down" chiral bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return The list of lines
     */
    std::vector<QLineF> calcDashedWedge(const QLineF& trimmed_line) const;

    /**
     * Calculate a wedge where lines zig-zag along the length of the wedge,
     * which is used to represent an "up or down" chiral bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return The list of lines for the zig-zags
     */
    std::vector<QLineF> calcSquigglyWedge(const QLineF& trimmed_line) const;

    /**
     * Calculate the lines needed to paint a double bond.  Note that bond
     * direction (i.e. a double "cis or trans" bond) is *not* taken into account
     * here.  For "cis or trans" double bonds (i.e. crossed double bonds),
     * calcSymmetricDoubleBondLines should be called directly instead of
     * calcDoubleBondLines.
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @param bond_line The untrimmed bond line
     * @return The two bond lines.  If the bond is aromatic, the second line of
     * the pair should be drawn dashed.
     */
    std::pair<QLineF, QLineF>
    calcDoubleBondLines(const QLineF& trimmed_line,
                        const QLineF& bond_line) const;

    /**
     * Determine whether this double bond should be drawn symmetrically (with
     * each line equally offset from the center line) or asymmetrically (with
     * one line along the center line and one line offset).  Note that bond
     * direction (i.e. a crossed "cis or trans" bond) is *not* taken into
     * account here.
     * @return true if the double bond is symmetrical, false otherwise
     */
    bool isDoubleBondSymmetric() const;

    /**
     * Calculate the lines needed to paint a symmetric double bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return The two bond lines
     */
    std::pair<QLineF, QLineF>
    calcSymmetricDoubleBondLines(const QLineF& trimmed_line) const;

    /**
     * Calculate the lines needed to paint an asymmetric double bond
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @param bond_line The untrimmed bond line
     * @return The two bond lines.  If the bond is aromatic, the second line of
     * the pair should be drawn dashed.
     */
    std::pair<QLineF, QLineF>
    calcAsymmetricDoubleBondLines(const QLineF& trimmed_line,
                                  const QLineF& bond_line) const;

    /**
     * Calculate the offset for an asymmetric double bond, i.e. where should
     * we draw the second bond line.
     * @param trimmed_line The bond line, trimmed so it doesn't overlap any atom
     * labels
     * @return A point giving the offset vector for the second bond line
     */
    QPointF calcAsymmetricDoubleBondOffset(const QLineF& trimmed_line) const;

    /**
     * For an asymmetric double bond, figure out which ring the second bond
     * should be drawn inside of.  If this bond is part of multiple rings, then
     * the criteria for picking the "best" one are:
     * - a ring with eight or fewer atoms is always preferred over a larger ring
     *   (i.e. avoid macrocycles if possible)
     * - after that, we want the ring with the most double bonds
     * - after that, we use the lowest ring index just to make the results
     *   reproducible
     * @param molecule The RDKit molecule
     * @param ring_info The RDKit ring info for the molecule
     * @return The ring index of the "best" ring.  If this bond is not part of
     * any rings, then -1 will be returned.
     */
    int findBestRingForDoubleBond(const RDKit::ROMol& molecule,
                                  const RDKit::RingInfo* ring_info) const;

    /**
     * Find the center of the specified ring
     * @param molecule The RDKit molecule
     * @param ring_info The RDKit ring info for the molecule
     * @param ring_index The index of the ring to find the center of
     * @return The ring center in local coordinates
     */
    QPointF findRingCenter(const RDKit::ROMol& molecule,
                           const RDKit::RingInfo* ring_info,
                           const int ring_index) const;

    /**
     * Determine whether the neighboring atoms of this bond (i.e. the neighbors
     * of the start and end atoms) are on the same side of the bond as offset,
     * or the opposite side.  In other words, if this bond was oriented
     * horizontally, would the neighboring atoms be above or below the bond?
     * @param molecule The RDKit molecule
     * @param neighbors_of The atom to check the neighbors of.  Should be either
     * the start or end atom of this bond.
     * @param atom_to_skip Skip this neighbor of neighbors_of.  Should be either
     * the end or start atom of this bond.
     * @param conf The 2d coordinates for atoms in the molecule
     * @param offset A point on one side of the bond or the other
     * @param line_endpoint A point along this bond.  The other end of the bond
     * is assumed to be at (0, 0) since we're working in a local coordinate
     * system centered about the start atom.
     * @param[in,out] num_same_side The count of neighbors on the same side of
     * the bond as offset.  This count will be updated for the neighbors of
     * neighbors_of.
     * @param[in,out] num_opposite_side The count of neighbors on the opposite
     * side of the bond as offset.  This count will be updated for the neighbors
     * of neighbors_of.
     */
    void countNeighboringAtomsBySide(
        const RDKit::ROMol& molecule, const RDKit::Atom* neighbors_of,
        const RDKit::Atom* atom_to_skip, const RDKit::Conformer& conf,
        const QPointF& offset, const QPointF& line_endpoint, int& num_same_side,
        int& num_opposite_side) const;

    /**
     * Determine whether two points are on the same or different sides of
     * a line.
     * @param point1 The first point
     * @param point2 The second point
     * @param line_endpoint A point on the line.  The other end of the line is
     * assumed to be (0, 0).
     * @return true if the points are on the same side of the line.  false
     * otherwise.
     */
    bool arePointsOnSameSideOfLine(const QPointF& point1, const QPointF& point2,
                                   const QPointF& line_endpoint) const;

    /**
     * Modify inner_line so that each end is at least 6 pixels shorter than
     * bond_line, unless the atom at that end of the bond is a terminal atom.
     * If inner_line has already been trimmed by 6 pixels or more (due to an
     * atom label), no additional trimming will be done on that end.
     * @param[in,out] inner_line The line to trim
     * @param offset  The offset used to calculate inner_line
     * @param bond_line The full length bond line (i.e. a line from the start
     * atom to the end atom)
     */
    void trimDoubleBondInnerLine(QLineF& inner_line, const QPointF& offset,
                                 const QLineF& bond_line) const;

    /**
     * Cross the given lines.  Used for "cis or trans" double bonds.
     * @param[in,out] line_one The first line.  One endpoint of this line will
     * be modified.
     * @param[in,out] line_two The second line.  One endpoint of this line will
     * be modified.
     */
    void crossDoubleBondLines(QLineF& line_one, QLineF& line_two) const;

    /**
     * Calculate two lines that are equally offset from a center line.
     * @param center_line The center line
     * @param offset The distance from the center line to each generated line
     * @return The pair of newly generated lines
     */
    std::pair<QLineF, QLineF> calcOffsetBondLines(const QLineF& center_line,
                                                  qreal offset) const;

    /**
     * Calculate the lines needed to paint a triple bond.
     *
     * @param trimmed_line A line from m_start_item to m_end_item that has been
     * trimmed so it doesn't overlap the atom item labels
     * @return The lines to paint.  This list will always contain exactly three
     * items.
     */
    std::vector<QLineF> calcTripleBondLines(const QLineF& trimmed_line) const;

    /**
     * Trim a line so that it doesn't overlap the atom item labels from either
     * bound atom.
     *
     * @param line The line to trim
     * @return The trimmed line
     */
    QLineF trimLineToBoundAtoms(const QLineF& line) const;

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
