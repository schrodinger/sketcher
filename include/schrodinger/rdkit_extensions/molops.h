/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: miscellaneous mol operations
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <boost/shared_ptr.hpp>
#include <vector>

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
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

/**
 * Adds all explicit hydrogens generating coordinates for each.
 * @param atom_ids vector of atom indexes where Hs should be added; if the
 *                  vector is empty, hydrogens will be added to the entire mol
 */
RDKIT_EXTENSIONS_API void addHs(RDKit::RWMol& rdk_mol,
                                std::vector<unsigned int> atom_ids = {});

/**
 * Removes hydrogens based on a common standard of what can be removed
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void removeHs(RDKit::RWMol& rdk_mol);

/**
 * Overload to remove only the Hs that are either included in the vector of
 * ints, or attached to atoms in the list.
 * @param rdk_mol rdkit mol
 * @param atom_ids vector of atom indexes where Hs should be removed
 *
 * @note atom_ids will be sorted and uniqued
 */
RDKIT_EXTENSIONS_API void removeHs(RDKit::RWMol& rdk_mol,
                                   std::vector<unsigned int> atom_ids);

//
// Helper api to extract a subgraph from an ROMol. Bonds, substance groups and
// stereo groups are only extracted to the subgraph if all participant atoms
// are selected by the `atom_ids` parameter.
//
// @param mol starting mol
// @param atom_ids the indices of atoms to extract. If an atom index falls
//                 outside of the acceptable atom indices, it is ignored.
// @param sanitize whether to sanitize the extracted mol.
//
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
ExtractMolFragment(const RDKit::ROMol& mol,
                   const std::vector<unsigned int>& atom_ids,
                   bool sanitize = true);

//
// Helper api to merge two monomeric mols. This api renames polymer ids in
// mol2 to prevent them from clashing with polymer ids in mol1.
//
// @param mol1 rdkit mol
// @param mol2 rdkit mol
//
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::ROMol>
CombineMonomericMols(const RDKit::ROMol& mol1, const RDKit::ROMol& mol2);

//
// Helper api to merge two mols with support for monomeristic mols.
//
// NOTE: If any of the inputs is a monomeristic mol, the output will be a
// monomeristic mol.
//
// @param mol1 rdkit mol
// @param mol2 rdkit mol
//
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::ROMol>
CombineMols(const RDKit::ROMol& mol1, const RDKit::ROMol& mol2);

} // namespace rdkit_extensions
} // namespace schrodinger
