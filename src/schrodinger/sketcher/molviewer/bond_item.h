#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include <QBrush>
#include <QGraphicsItem>
#include <QLineF>
#include <QPen>
#include <QRectF>
#include <QtGlobal>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

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
 * A Qt graphics item for representing atomistic bonds in a molviewer Scene.
 */
class SKETCHER_API BondItem : public AbstractBondOrConnectorItem
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
    BondItem(const RDKit::Bond* bond, const AtomItem& start_atom_item,
             const AtomItem& end_atom_item, const Fonts& fonts,
             const BondDisplaySettings& settings,
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

    /**
     * @return the bond associated with this item
     */
    const RDKit::Bond* getBond() const;

  protected:
    const AtomItem& m_start_item;
    const AtomItem& m_end_item;
    std::vector<ToPaint> m_to_paint;
    QPen m_solid_pen;
    QPen m_dashed_pen;
    QPen m_chirality_pen;
    QBrush m_solid_brush = QBrush(Qt::black);
    qreal m_text_angle;
    QPointF m_text_pos;
    QSizeF m_text_size;

    const Fonts& m_fonts;

    const BondDisplaySettings& m_settings;

    QString m_annotation_text;

    std::vector<QColor> m_colors;

    /**
     * Get tooltip text based on bond stereochemistry or query labels.
     * Returns tooltip text for features related to stereochemistry or query
     * features. If a bond has both stereochemistry and query features, both are
     * shown on separate lines (e.g., "Stereo: (E)\nQuery: Any bond type").
     * @return tooltip text, or empty string if no tooltip should be shown
     */
    QString getTooltip() const;

    /**
     * @return a list of the colors to be used for painting this bond. If the
     * user has selected a color for this bond, then the list will contain only
     * that. If the user has not selected a color for the bond, but has selected
     * colors for the atoms, then the list will contain both (if they differ)
     * and the bond color will change at its half point. If the user has not
     * selected a color, m_setting.color will be used.
     */
    std::vector<QColor> getColors() const;

    /**
     * paint all bond lines and polygons using the given painter. This function
     * is called by paint() before painting the bond annotation. If an
     * annotation is present the bond lines need to be partially transparent
     * behind it to make the annotation readable, so this function will be
     * called twice with different painter opacities and using the annotation
     * polygon as a clipping region.
     */
    void paintBondLinesAndPolygons(QPainter* painter);

    /**
     * Calculate the lines and polygons needed to paint this bond.  Note that
     * this method does *not* do any actual painting.  (The output of this
     * method should be stored in m_to_paint, and the actual painting will be
     * done in `paint`.)
     *
     * @param bond_line A line from m_start_item to m_end_item
     * @param bond_type The type of bond to paint
     * @return All of the lines and polygons to paint, along with the pens to be
     * used for painting.
     */
    std::vector<ToPaint>
    calculateLinesToPaint(const QLineF& bond_line,
                          const RDKit::Bond::BondType bond_type) const;

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
     * For an asymmetric  bond (double or aromatic), figure out which ring the
     * second bond should be drawn inside of.  If this bond is part of multiple
     * rings, then the criteria for picking the "best" one are:
     * - a ring with eight or fewer atoms is always preferred over a larger ring
     *   (i.e. avoid macrocycles if possible)
     * - after that, we want the ring with the most double or aromatic bonds
     * - after that, we use the lowest ring index just to make the results
     *   reproducible
     * @param molecule The RDKit molecule
     * @param ring_info The RDKit ring info for the molecule
     * @return The ring index of the "best" ring.  If this bond is not part of
     * any rings, then -1 will be returned.
     */
    int findBestRingForBond(const RDKit::ROMol& molecule,
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
     * Determine the parameters to be used for painting the stereo or query
     * annotation.  The label is always drawn parallel to the bond.
     * @param label The label text to be displayed.
     * @param draw_text_above_bond Whether to draw the text above the bond
     * or over the bond
     * @return The calculated parameters:
     * - angle: the angle of the annotation text
     * - text_pos: the position of center of the label,
     * - text_size: size of the bounding rect of the label
     */
    std::tuple<qreal, QPointF, QSizeF>
    getStereoAnnotationParameters(const QString& label,
                                  const bool draw_text_above_bond) const;

    /** return the annotation bounding polygon for the bond. Note that the text
     * could be rotated, so we need a polygon rather than a rectangle. The text
     * painted by paintAnnotation() is guaranteed to stay within this polygon.
     */
    QPolygonF getAnnotationPolygon();

    /** Paint text annontation .
     *
     * @param painter  the qt painter to be used
     * @param angle  the angle of the annotation text
     * @param text_pos the position of center of the label,
     * @param text_size size of the label
     * @param text  the text to be displayed
     */
    void paintAnnotation(QPainter* painter, qreal angle,
                         const QPointF& text_pos, const QSizeF& text_size,
                         const QString& text);

    /** @return the stereo label associated with this bond
     */
    QString getStereoLabel();
};

} // namespace sketcher
} // namespace schrodinger
