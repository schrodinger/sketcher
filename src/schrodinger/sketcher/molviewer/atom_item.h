/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <vector>

#include <QPen>
#include <QRectF>
#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

class QGraphicsItem;

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{
namespace sketcher
{

enum class HsDirection { RIGHT, LEFT, UP, DOWN };

/**
 * A Qt graphics item for representing atoms in a molviewer Scene.  This class
 * is responsible for painting all aspects of the atom label, which can include
 *
 * - the main text label for the atom.  This normally contains the atom's
 *   element, but can also contain the R group or a query string.
 * - implicit hydrogens
 * - lone pair indicator
 * - valence error indicator
 *
 * as well as a handful of additional aspects, depending on the atom.  Note that
 * labels are frequently not shown for certain carbons, so this class won't
 * paint anything in those scenarios.  Also note that this class is *not*
 * responsible for drawing bonds.  See BondItem for that.
 */
class SKETCHER_API AtomItem : public AbstractGraphicsItem
{
  public:
    /**
     * Note that this class does not take ownership of atom, fonts, or settings.
     * These objects must not be destroyed while the AtomItem is in use.
     *
     * @param atom The RDKit atom that this item should represent.
     *
     * @param fonts The fonts to be used when painting this item.  To change the
     * fonts, the calling code (i.e. the Scene) must update this object and then
     * call updateCachedData.
     *
     * @param settings The settings to be used for painting this item.  To
     * change settings, the calling code (i.e. the Scene) must update this
     * object and then call updateCachedData.
     *
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     * @pre atom != nullptr
     * @pre atom->hasOwningMol()
     */
    AtomItem(RDKit::Atom* atom, Fonts& fonts, AtomItemSettings& settings,
             QGraphicsItem* parent = nullptr);

    // Type and type() are required for qgraphicsitem_cast.  Note that this enum
    // implementation is based on the sample code from the QGraphicsItem::Type
    // documentation.
    enum { Type = static_cast<int>(ItemType::ATOM) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

    /**
     * Return a list of bounding rectangles for all individual components of
     * this item.  These subrects are used to precisely trim the bond lines so
     * that they don't overlap with the atom labels.  Note that this method
     * returns subrects starting from the center of the atom and working
     * outwards.  This behavior is required for proper functioning of bond line
     * trimming due to assumptions made in BondItem::trimLineToRect.
     */
    const std::vector<QRectF>& getSubrects() const;

    /**
     * set all label rects to invalid rects and all label text to empty strings
     */
    void clearLabels();

    /**
     * create a label for the isotope number (if present)
     */
    void updateIsotopeLabel();

    /**
     * create a label for the charge and radical information  (if present)
     */
    void updateChargeAndRadicalLabel();

    /**
     * create a label for the hydrogens on this atom (if any)
     */
    void updateHsLabel();

    /**
     * position all the labels around the atom. Hs can be to the right, left,
     * top or bottom depending on the binding pattern
     */
    void positionLabels();

    /**
     * create and position labels (dots) for the unpaired electrons on this atom
     * (if any)
     */
    void updateUnpairedElectronsLabels();

    /**
     * Return whether the label for this atom is visible.  (Labels for some
     * carbons may be hidden depending on the atom item settings.)
     */
    bool labelIsVisible() const;

    /**
     * @param avoid_subrects whether or not to consider the atom's subrects
     * in the calculation: when picking a position for something close to the
     * atom (e.g. chiral labels), this make sure that the returned position
     * doesn't clash with the atom labels. When picking positions for something
     * far away from the label (e.g. another atom bound to this), this should be
     * false
     * @return a position on the canvas (in coordinates relative to this item)
     * that is near this atom but doesn't clash with it or its neighbors.
     */
    QPointF findPositionInEmptySpace(bool avoid_subrects) const;

  protected:
    /**
     * @Return the direction Hs labels would be drawn on this atom if the label
     * is shown and it has hydrogens
     */
    HsDirection findHsDirection() const;

    // Creating a shared_ptr to an RDKit Atom (or Bond) implicitly creates a
    // copy of the Atom, which means that the new Atom is no longer part of the
    // original molecule, which leads to problems.  Because of this, we store a
    // raw pointer instead.  The RDKit molecule takes care of the lifetime of
    // the Atom, so an AtomItem instance must be deleted as soon as its
    // associated Atom is deleted.
    //
    // Also note that m_atom should only be accessed from within
    // updateCachedData to ensure that we can properly notify the scene of any
    // AtomItem changes *before* they happen.
    RDKit::Atom* const m_atom;
    QString m_main_label_text;
    QRectF m_main_label_rect;
    QRectF m_H_label_rect;
    QRectF m_H_count_label_rect;
    QString m_H_count_label_text;
    QRectF m_charge_and_radical_rect;
    QString m_charge_and_radical_label_text;
    QString m_isotope_label_text;
    QRectF m_isotope_rect;

    std::vector<QRectF> m_subrects;
    Fonts& m_fonts;
    AtomItemSettings& m_settings;
    bool m_label_is_visible;
    QPen m_pen;
    QPen m_valence_error_pen;
    QBrush m_valence_error_brush;

    bool determineLabelIsVisible() const;

    /**
     * @Return whether the associated atom has a valence error that should be
     * displayed
     */
    bool shouldDisplayValenceError() const;
};

/**
 * @param xyz RDKit atom coordinates to transform
 * @return resulting 2D scene coordinates
 *
 * TODO: Move this to another header/source during SKETCH-1841
 */
SKETCHER_API QPointF to_scene_xy(const RDGeom::Point3D& xyz);

} // namespace sketcher
} // namespace schrodinger
