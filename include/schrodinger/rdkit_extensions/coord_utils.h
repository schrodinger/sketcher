#pragma once

#include <vector>

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

/**
 * Adds a 2D conformer if none present using the above function
 */
RDKIT_EXTENSIONS_API void update_2d_coordinates(RDKit::ROMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
