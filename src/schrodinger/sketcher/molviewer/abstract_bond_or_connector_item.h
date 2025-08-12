#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

class QPointF;

namespace RDKit
{
class Bond;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * An abstract parent class for all graphics items that represent bonds or
 * monomer connectors (which are stored as RDKit::Bonds).
 */
class SKETCHER_API AbstractBondOrConnectorItem : public AbstractGraphicsItem
{
  public:
    AbstractBondOrConnectorItem(const RDKit::Bond* bond,
                                QGraphicsItem* parent = nullptr);

    /**
     * @return the RDKit bond associated with this item
     */
    const RDKit::Bond* getBond() const;

    /**
     * @return the mid point of the bond (in local coordinates)
     */
    QPointF getMidpoint() const;

  protected:
    // Creating a shared_ptr to an RDKit Bond (or Atom) implicitly creates a
    // copy of the Bond, which means that the new Bond is no longer part of the
    // original molecule, which leads to problems. Because of this, we store a
    // raw pointer instead. The RDKit molecule takes care of the lifetime of the
    // Bond, so a BondItem instance must be deleted as soon as its associated
    // Bond is deleted.
    //
    // Also note that m_bond should only be accessed from within
    // updateCachedData to ensure that we can properly notify the scene of any
    // BondItem changes *before* they happen.
    const RDKit::Bond* const m_bond;

    /**
     * The midpoint of the bond, which is used to determine whether the bond is
     * included in a marquee selection.  Subclasses are responsible for setting
     * this value in updateCachedData.
     */
    QPointF m_midpoint;
};

} // namespace sketcher
} // namespace schrodinger
