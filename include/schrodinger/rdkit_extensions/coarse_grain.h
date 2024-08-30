/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: tools for Coarse-grain ROMols
 *
 * A "coarse grain" (CG) ROMol uses RDKit atoms to represent monomers. Chains
 * are represented by COP Substance Groups on the ROMol.
 *
 * For use with functionality in schrodinger::rdkit_extensions
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <string_view>

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
class ROMol;
class RWMol;
class Atom;
class Bond;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger::rdkit_extensions
{

using Monomer = RDKit::Atom;

enum class ChainType { PEPTIDE, RNA, DNA, CHEM };
enum class ConnectionType { FORWARD, SIDECHAIN };
enum class MonomerType { REGULAR, /* LIST, WILDCARD, */ SMILES };

RDKIT_EXTENSIONS_API ChainType to_chain_type(std::string_view chain_type);

RDKIT_EXTENSIONS_API std::string to_string(ChainType chain_type);

/*
 * Add a monomer to the molecule
 *
 * @param cg_mol The CG to add the monomer to
 * @param name The name of the monomer
 * @param residue_number The residue number of the monomer
 * @param chain_id The chain ID of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_EXTENSIONS_API size_t add_monomer(
    RDKit::RWMol& cg_mol, std::string_view name, int residue_number,
    std::string_view chain_id, MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a monomer to the molecule. Overload that uses the last monomer
 * added to the molecule to determine the chain ID and residue number.
 *
 * @param cg_mol The CG to add the monomer to
 * @param name The name of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_EXTENSIONS_API size_t
add_monomer(RDKit::RWMol& cg_mol, std::string_view name,
            MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a connection between two monomers in the molecule. The connection has
 * directionality that starts at monomer1 and ends at monomer2.
 *
 * @param mol The molecule to add the connection to
 * @param monomer1 The index of the first monomer
 * @param monomer2 The index of the second monomer
 * @param connection_type The type of connection to add
 */
RDKIT_EXTENSIONS_API void
add_connection(RDKit::RWMol& mol, size_t monomer1, size_t monomer2,
               ConnectionType connection_type = ConnectionType::FORWARD);

// overload for helm writer
RDKIT_EXTENSIONS_API void add_connection(RDKit::RWMol& mol, size_t monomer1,
                                         size_t monomer2,
                                         const std::string& linkage);

// Discards existing chains and reassigns monomers to sequential chains.
// (in HELM world, "chains" are called "polymers")
RDKIT_EXTENSIONS_API void assign_chains(RDKit::RWMol& mol);

} // namespace schrodinger::rdkit_extensions