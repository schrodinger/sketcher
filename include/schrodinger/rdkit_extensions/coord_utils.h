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
 * If no 2d conformer is present on the given molecule, add one using
 * compute2DCoords above.  If a 2d conformer is present, ensure that bonds are
 * RDDepict::BOND_LEN units long, rescaling if necessary.
 */
RDKIT_EXTENSIONS_API void update_2d_coordinates(RDKit::ROMol& mol);

/**
 * Determine the length of a typical bond in the given molecule.  If the
 * molecule has no coordinates or no bonds, -1 will be returned.
 */
RDKIT_EXTENSIONS_API double
get_most_common_bond_length(const RDKit::ROMol& mol);

/**
 * Ensure that bonds in the given molecule are RDDepict::BOND_LEN units long and
 * rescale if necessary.  Nothing will be done if the molecule has no
 * coordinates or no bonds.
 */
RDKIT_EXTENSIONS_API void rescale_bond_length_if_needed(RDKit::ROMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
