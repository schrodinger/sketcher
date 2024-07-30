/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: atomistic ROMol -> CG mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <boost/shared_ptr.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

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
 * Converts an atomistic ROMol into a CG ROMol that can
 * be written to HELM format using the HELM writer.
 *
 * @param atomistic_mol Atomistic molecule to convert to CG
 * @param use_residue_info Whether to use PDBAtomResidueInfo to determine
 * monomer boundaries
 * @return CG molecule
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
atomistic_to_cg(const RDKit::ROMol& atomistic_mol,
                bool use_residue_info = false);

/**
 * Identify monomers within an atomistic molecule
 *
 * For testing purposes.
 */
RDKIT_EXTENSIONS_API std::vector<std::vector<int>>
get_monomers(const RDKit::ROMol& mol);

/**
 * Build an atomistic molecule from a CG molecule
 *
 * For now, assumes that the CG molecule was built using the HELM parser.
 *
 * @param cg_mol CG molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
cg_to_atomistic(const RDKit::ROMol& cg_mol);

} // namespace rdkit_extensions
} // namespace schrodinger
