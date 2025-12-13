#pragma once

#include <unordered_set>
#include <vector>

#include "schrodinger/sketcher/definitions.h"

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
namespace sketcher
{

/**
 * Ensure that the specified molecule contains exactly one 2d conformer. If no
 * 2d conformers are present, one will be added using compute2DCoords. If
 * a 2d conformer is present, bonds will be rescaled as needed to ensure that
 * they are RDDepict::BOND_LEN units long. If the input structure contains a 3d
 * conformer, it will be kept in addition to the 2d conformer. Any additional 2d
 * or 3d conformers will be discarded.
 */
SKETCHER_API void update_2d_coordinates(RDKit::ROMol& mol);

/**
 * Determine the length of a typical bond in the given molecule.  If the
 * molecule has no coordinates or no bonds, -1 will be returned.
 */
SKETCHER_API double get_most_common_bond_length(const RDKit::ROMol& mol);

/**
 * Ensure that bonds in the given molecule are RDDepict::BOND_LEN units long and
 * rescale if necessary.  Nothing will be done if the molecule has no
 * coordinates or no bonds.
 */
SKETCHER_API void rescale_bond_length_if_needed(RDKit::ROMol& mol);

} // namespace sketcher
} // namespace schrodinger
