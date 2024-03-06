#pragma once

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class ROMol;
class Conformer;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Fix the coordinates for any variable attachment bonds in the given molecule.
 * This should be run after calling rdkit_extensions::compute2DCoords on any
 * molecule with variable attachment bonds, since RDKit's minimzer doesn't
 * understand variable attachment bonds and thinks they're not bound to the
 * variable attachment atoms.
 */
SKETCHER_API void fix_variable_attachment_bond_coordinates(RDKit::ROMol& mol);

} // namespace sketcher
} // namespace schrodinger
