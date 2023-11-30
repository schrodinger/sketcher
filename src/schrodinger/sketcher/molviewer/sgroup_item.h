#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include <GraphMol/SubstanceGroup.h>

class QPainter;
namespace RDKit
{
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics item for painting substance groups in a molviewer Scene.
 */
class SKETCHER_API SGroupItem : public AbstractGraphicsItem
{

  public:
    /**
     * Note that this class does not take ownership of the
     * atom items.  These objects must not be destroyed while the SGroupItem is
     * in use.
     *
     * @param sgroup The RDKit substance group that this item should represent.
     * @param fonts to render the labels
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     */
    SGroupItem(const RDKit::SubstanceGroup& sgroup, const Fonts& fonts,
               QGraphicsItem* parent = nullptr);

    const RDKit::SubstanceGroup* getSubstanceGroup() const;

    // reimplemented methods
    // Type and type() are required for qgraphicsitem_cast.  Note that this
    // enum implementation is based on the sample code from the
    // QGraphicsItem::Type documentation.
    enum { Type = static_cast<int>(ItemType::S_GROUP) };
    int type() const override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
    void updateCachedData() override;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;

  private:
    /**
     * @return the vector that is perpendicular to the line between point1 and
     * point2 (the two endpoints of a bracket's long side) and points towards
     * the inside of the sgroup
     */
    RDGeom::Point3D
    computeShortSideForBracket(const RDGeom::Point3D& point1,
                               const RDGeom::Point3D& point2) const;

    /**
     * @return the predicted highlighting path for this item. This highlights
     * the brackets, the labels, and the atoms and bonds in the sgroup.
     */

    QPainterPath computePredictiveHighlightingPath() const;

    const RDKit::SubstanceGroup& m_sgroup;
    QPainterPath m_brackets_path;
    const Fonts& m_fonts;

    /**
     * @return data required to position the bracket labels.  This data consists
     * of a QLineF and a QPointF in Scene coordinates:
     * - the long side of the rightmost bracket.  Both labels are positioned
     *   relative to this bracket.
     * - the offset vector for the labels.  The labels are positioned by adding
     *   this offset to either end of the rightmost bracket.
     */
    std::pair<QLineF, QPointF> getPositionsForLabels() const;

    /**
     * @return The QPainterPath for the brackets around the substance group.
     * Each bracket goes through the middle of a connection bond and is
     * perpendicular to it. A connection bond is any bond that connects an
     * atom inside the substance group to an atom outside of it.
     */
    QPainterPath getBracketPath() const;

    QString m_label;
    QString m_repeat;

    QRectF m_untransformed_label_rect;
    QRectF m_untransformed_repeat_rect;
    QTransform m_labels_transform;
};

} // namespace sketcher
} // namespace schrodinger
