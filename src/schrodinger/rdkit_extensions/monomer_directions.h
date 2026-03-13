#pragma once

#include <rdkit/Geometry/point.h>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger::rdkit_extensions
{

enum class Direction { N, S, E, W, NW, NE, SW, SE };

/**
 * Convert a Direction to a unit vector in RDKit mol coordinates.
 * In mol coordinates, +X is right and +Y is up.
 */
RDKIT_EXTENSIONS_API RDGeom::Point3D direction_to_unit_vector(Direction dir);

} // namespace schrodinger::rdkit_extensions
