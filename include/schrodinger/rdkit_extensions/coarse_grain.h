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
using Chain = RDKit::SubstanceGroup;

enum class ChainType { PEPTIDE, RNA, DNA, CHEM };
enum class ConnectionType { FORWARD, SIDECHAIN };
enum class MonomerType { REGULAR, /* LIST, WILDCARD, */ SMILES };

/*
 * Add a monomer to the molecule
 *
 * @param chain The chain to add the monomer to
 * @param name The name of the monomer
 * @param monomer_type The type of monomer to add
 */
RDKIT_EXTENSIONS_API size_t
add_monomer(Chain& chain, std::string_view name,
            MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a connection between two monomers in the molecule
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

/*
 * Create a new chain (COP substance group) and add it to the molecule
 */
RDKIT_EXTENSIONS_API Chain& add_chain(RDKit::RWMol&, ChainType chain_type);

// overload for helm writer
RDKIT_EXTENSIONS_API Chain& add_chain(RDKit::RWMol&,
                                      std::string_view chain_name);

// Discards existing chains and reassigns monomers to sequential chains.
// (in HELM world, "chains" are called "polymers")
RDKIT_EXTENSIONS_API void assign_chains(RDKit::RWMol& mol);

} // namespace schrodinger::rdkit_extensions