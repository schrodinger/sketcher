#pragma once

#include <rdkit/Geometry/point.h>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger::rdkit_extensions
{

enum class Direction { N, S, E, W, NW, NE, SW, SE };

/**
 * Convert a Direction to a vector in RDKit mol coordinates. Vectors in cardinal
 * directions will be 1 unit long, and vectors in diagonal directions will have
 * positive or negative 1 for their X and Y coordinates (i.e. they'll be
 * slightly longer that the cardinal vectors).  Note that in RDKit coordinates,
 * +X is right and +Y is up (unlike Qt coordinates).
 */
RDKIT_EXTENSIONS_API RDGeom::Point3D direction_to_vector(Direction dir);

} // namespace schrodinger::rdkit_extensions
