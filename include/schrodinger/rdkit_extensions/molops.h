/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: miscellaneous mol operations
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
class Atom;
class ROMol;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Level of RDKit sanitization to apply
 *
 * FULL indicates all possible sanitization routines
 *
 * PARTIAL indicates only the following MolOps:: operations are performed:
 * - clear out any cached properties on atoms/bonds and update computed ones
 * - symmetrize SSSR
 * - set conjugation
 * - set hybridization
 * - adjust Hydrogen counts
 */
enum class Sanitization { FULL, PARTIAL };

/**
 * @param mol rdkit mol
 * @param sanitization what level of sanitization to apply
 */
RDKIT_EXTENSIONS_API void
apply_sanitization(RDKit::RWMol& mol,
                   Sanitization sanitization = Sanitization::FULL);

/*
 * Assigns stereochemistry to the given molecule, from 3D or 2D if possible
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void assign_stereochemistry(RDKit::ROMol& mol);

/**
 * @param atom rdkit atom
 * @return whether the atom is an attachment point
 */
RDKIT_EXTENSIONS_API bool is_attachment_point_dummy(const RDKit::Atom& atom);

/**
 * Set enhanced stereo for any chiral atoms that don't already have it.
 * If the chiral flag is on, ungrouped chiral atoms will go into the ABS
 * group. If the flag is off or not present, ungrouped atoms will go into
 * a new AND group.
 */
RDKIT_EXTENSIONS_API void
add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol);

/**
 * Set bond directions based on the input data from SDF
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void reapply_molblock_wedging(RDKit::ROMol& rdk_mol);

/**
 * Removes hydrogens based on a common standard of what can be removed
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void removeHs(RDKit::RWMol& rdk_mol);

} // namespace rdkit_extensions
} // namespace schrodinger
