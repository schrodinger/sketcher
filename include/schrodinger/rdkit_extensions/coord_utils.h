#pragma once

#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace RDKit
{
class Atom;
class Conformer;
class ROMol;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @return the centroid of a set of points
 *
 */
RDKIT_EXTENSIONS_API RDGeom::Point3D
compute_centroid(const std::vector<RDGeom::Point3D>& positions);

/**
 * Generate 2D coordinates for the given molecule using RDDepict with ring
 * templating enabled; note this function forces RDKit coordinate generation
 * and will ignore SetPreferCoordGen if set.
 *
 * @param mol rdkit mol
 * @param frozen_ids vector of atom indexes to NOT generate coordinates for
 * @return ID of the conformation added to the molecule containing the 2D coords
 */
RDKIT_EXTENSIONS_API unsigned int
compute2DCoords(RDKit::ROMol& mol,
                const std::vector<unsigned int>& frozen_ids = {});

} // namespace rdkit_extensions
} // namespace schrodinger
