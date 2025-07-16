#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

namespace RDKit
{
class Atom;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * An abstract parent class for all graphics items that represent atoms or
 * monomers (which are stored as RDKit::Atoms).
 */
class SKETCHER_API AbstractAtomOrMonomerItem : public AbstractGraphicsItem
{
  public:
    AbstractAtomOrMonomerItem(const RDKit::Atom* atom_or_monomer,
                              QGraphicsItem* parent = nullptr);

    /**
     * @return the RDKit atom associated with this item
     */
    const RDKit::Atom* getAtom() const;

  protected:
    // Creating a shared_ptr to an RDKit Atom (or Bond) implicitly creates a
    // copy of the Atom, which means that the new Atom is no longer part of the
    // original molecule, which leads to problems.  Because of this, we store a
    // raw pointer instead.  The RDKit molecule takes care of the lifetime of
    // the Atom, so the graphics item instance must be deleted as soon as its
    // associated Atom is deleted.
    //
    // Also note that m_atom should only be accessed from within
    // updateCachedData to ensure that we can properly notify the scene of any
    // changes *before* they happen.

    const RDKit::Atom* const m_atom;
};

} // namespace sketcher
} // namespace schrodinger
