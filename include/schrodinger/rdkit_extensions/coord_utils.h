#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

namespace RDKit
{
class Conformer;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Generate 2D coordinates for the given molecule using RDDepict with ring
 * templating enabled
 *
 * @return ID of the conformation added to the molecule containing the 2D coords
 */
RDKIT_EXTENSIONS_API unsigned int compute2DCoords(RDKit::ROMol& mol);

/**
 * Adds a 2D conformer if none present
 */
RDKIT_EXTENSIONS_API void update_2d_coordinates(RDKit::ROMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
