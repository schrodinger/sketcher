/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: atomistic ROMol -> Monomer Mol
 conversion
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
 * Converts an atomistic ROMol into a monomeric ROMol that can
 * be written to HELM format using the HELM writer.
 *
 * @param atomistic_mol Atomistic molecule to convert to monomeric
 * @param try_residue_info Whether to first try using PDBAtomResidueInfo to
 * determine monomer boundaries. If set to false, will always use SMARTS
 * matching method to identify monomers.
 * @return monomeric molecule
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
toMonomeric(const RDKit::ROMol& atomistic_mol, bool try_residue_info = true);

/**
 * Identify monomers within an atomistic molecule
 *
 * For testing purposes.
 */
RDKIT_EXTENSIONS_API std::vector<std::vector<unsigned int>>
getMonomers(const RDKit::ROMol& mol);

/**
 * Build an atomistic molecule from a monomeric molecule
 *
 * For now, assumes that the monomeric molecule was built using the HELM
 * parser.
 *
 * @param monomer_mol monomeric molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
toAtomistic(const RDKit::ROMol& monomer_mol);

} // namespace rdkit_extensions
} // namespace schrodinger
