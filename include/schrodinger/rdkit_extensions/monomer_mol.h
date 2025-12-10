/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: tools for Monomeric ROMols
 *
 * A "MonomerMol" is an ROMol that uses RDKit atoms to represent monomers.
 Chains
 * are represented via the PDBAtomResidueInfo structs on atoms, and linkages are
 * stored as a LINKAGE property on bonds in the form of RX-RY, where X is the
 attachment
 * point used on the begin monomer and Y is the attachment point used on the
 monomer.
 *
 * For use with functionality in schrodinger::rdkit_extensions
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>
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

RDKIT_EXTENSIONS_API ChainType toChainType(std::string_view chain_type);

RDKIT_EXTENSIONS_API std::string toString(ChainType chain_type);

/*
 * Add a monomer to the molecule
 *
 * @param monomer_mol The monomeric molecule to add the monomer to
 * @param name The name of the monomer
 * @param residue_number The residue number of the monomer
 * @param chain_id The chain ID of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_EXTENSIONS_API size_t addMonomer(
    RDKit::RWMol& monomer_mol, std::string_view name, int residue_number,
    std::string_view chain_id, MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a monomer to the molecule. Overload that uses the last monomer
 * added to the molecule to determine the chain ID and residue number.
 *
 * @param monomer_mol The monomeric molecule to add the monomer to
 * @param name The name of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_EXTENSIONS_API size_t
addMonomer(RDKit::RWMol& monomer_mol, std::string_view name,
           MonomerType monomer_type = MonomerType::REGULAR);

/**
 * Create and return a new monomer
 * @param name The name of the monomer
 * @param chain_id The chain ID of the monomer
 * @param residue_number The residue number of the monomer
 * @param is_smiles Whether the monomer is a SMILES monomer
 */
RDKIT_EXTENSIONS_API std::unique_ptr<Monomer>
makeMonomer(const std::string_view name, const std::string_view chain_id,
            const int residue_number, const bool is_smiles);

/*
 * Mutate a monomer in the provided monomer mol to a new monomer specified
 * by a HELM symbol. Assumes the monomer is mutated to another monomer of the
 * same chain type (such as PEPTIDE or RNA).
 *
 * @param monomer_mol The monomeric molecule to mutate
 * @param monomer_idx The index of the monomer to mutate
 * @param helm_symbol The HELM symbol to mutate the monomer to
 */
RDKIT_EXTENSIONS_API void mutateMonomer(RDKit::ROMol& monomer_mol,
                                        unsigned int monomer_idx,
                                        std::string_view helm_symbol);

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
addConnection(RDKit::RWMol& mol, size_t monomer1, size_t monomer2,
              ConnectionType connection_type = ConnectionType::FORWARD);

// overload for helm writer
RDKIT_EXTENSIONS_API void addConnection(RDKit::RWMol& mol, size_t monomer1,
                                        size_t monomer2,
                                        const std::string& linkage,
                                        const bool is_custom_bond = false);

// Discards existing chains and reassigns monomers to sequential chains.
// (in HELM world, "chains" are called "polymers")
RDKIT_EXTENSIONS_API void assignChains(RDKit::RWMol& mol);

} // namespace schrodinger::rdkit_extensions
