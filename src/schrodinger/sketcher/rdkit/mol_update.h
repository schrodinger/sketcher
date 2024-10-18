#pragma once

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class ROMol;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Prepare an RDKit molecule for use in the sketcher, including generation
 * of 2D coordinates and certain updates dealing with stereochemistry.
 *
 * @param mol The molecule to prepare for use in the sketcher
 */
SKETCHER_API void prepare_mol(RDKit::ROMol& mol);

/**
 * Update an RDKit molecule after any change is made to its underlying
 * chemistry. This includes any add/update/delete operations on atoms, bonds,
 * and sgroups.
 *
 * @param mol The molecule to update
 */
SKETCHER_API void update_molecule_on_change(RDKit::RWMol& mol);

} // namespace sketcher
} // namespace schrodinger
