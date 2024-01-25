#include "schrodinger/sketcher/model/non_molecular_object.h"

namespace schrodinger
{
namespace sketcher
{

NonMolecularObject::NonMolecularObject(const NonMolecularType& type,
                                       const RDGeom::Point3D& coords,
                                       const NonMolecularTag tag) :
    m_type(type),
    m_coords(coords),
    m_tag(tag)
{
}

NonMolecularType NonMolecularObject::getType() const
{
    return m_type;
}

RDGeom::Point3D NonMolecularObject::getCoords() const
{
    return m_coords;
}

void NonMolecularObject::setCoords(const RDGeom::Point3D& coords)
{
    m_coords = coords;
}

NonMolecularTag NonMolecularObject::getTag() const
{
    return m_tag;
}

} // namespace sketcher
} // namespace schrodinger
