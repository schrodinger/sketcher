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

/**
 * Identify the residues in a molecule and set the AtomMonomerInfo for each
 * atom, including atom name, residue name and number, and chain ID.
 *
 * The atom names come from the CXSMILES data in the SMILES field of the
 * monomer database; the residue names come from the PDBCODE field. Chain names
 * are assigned as A, B, C, etc. in an unspecified order, but only for chains
 * containing more than one monomer; single-monomer chains are assigned a space
 * character as the chain ID. Residue numbers start from 1 for each chain.
 */
RDKIT_EXTENSIONS_API void addResidueInfo(RDKit::ROMol& atomistic_mol);

} // namespace rdkit_extensions
} // namespace schrodinger
