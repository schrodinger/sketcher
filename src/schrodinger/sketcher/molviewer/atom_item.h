/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <vector>

#include <QPen>
#include <QRectF>
#include <QString>

#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

class QGraphicsItem;

namespace RDKit
{
class Atom;
}

namespace schrodinger
{
namespace sketcher
{

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
class AtomItem : public AbstractGraphicsItem
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

  protected:
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
    std::vector<QRectF> m_subrects;
    Fonts& m_fonts;
    AtomItemSettings& m_settings;
    bool m_label_is_visible;
    QPen m_pen;

    bool determineLabelIsVisible() const;
};

} // namespace sketcher
} // namespace schrodinger
