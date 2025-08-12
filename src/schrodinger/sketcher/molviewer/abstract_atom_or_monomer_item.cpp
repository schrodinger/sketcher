#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"

namespace schrodinger
{
namespace sketcher
{

AbstractAtomOrMonomerItem::AbstractAtomOrMonomerItem(
    const RDKit::Atom* atom_or_monomer, QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_atom(atom_or_monomer)
{
}

const RDKit::Atom* AbstractAtomOrMonomerItem::getAtom() const
{
    return m_atom;
}

} // namespace sketcher
} // namespace schrodinger
