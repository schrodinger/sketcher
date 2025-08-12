#pragma once

#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/tags.h"

namespace schrodinger
{
namespace sketcher
{

class MolModel;

/**
 * The accepted types of non-molecular objects
 */
enum class NonMolecularType {
    RXN_ARROW,
    RXN_PLUS,
};

/**
 * A class to represent a non-molecular object: an object that is (undoably)
 * stored in MolModel, but that isn't part of the RDKit molecule.  This includes
 * reaction arrows and plus signs.
 */
class SKETCHER_API NonMolecularObject
{
  public:
    /**
     * @param type The type of non-molecular object
     * @param coords Coordinates for the center of the object.  Note that these
     * coordinates use the same coordinate system as the RDKit molecule.
     * @param tag The unique tag used to identify this non-molecular object
     */
    NonMolecularObject(const NonMolecularType& type,
                       const RDGeom::Point3D& coords,
                       const NonMolecularTag tag);

    NonMolecularType getType() const;

    RDGeom::Point3D getCoords() const;

  protected:
    NonMolecularType m_type;
    RDGeom::Point3D m_coords;
    NonMolecularTag m_tag;

    NonMolecularTag getTag() const;

    void setCoords(const RDGeom::Point3D& coords);

    friend class MolModel; // MolModel needs access to getTag and setCoords
};

} // namespace sketcher
} // namespace schrodinger
