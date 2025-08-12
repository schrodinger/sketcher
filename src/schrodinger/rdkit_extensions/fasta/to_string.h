#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <string>

namespace RDKit
{
class ROMol;
}

namespace fasta
{

/*
 * Converts a monomeric mol to a fasta input. The monomeric mol can
 * contain multiple chains and the chains can be heterogenous.
 *
 * @param mol: the monomeric mol
 * @throws std::invalid_argument: If the input is malformed
 *
 *
 * NOTE: Unsupported features:
 *      * atomistic mols
 *      * SRU substance groups
 *      * non-peptide and non-nucleotide polymers
 *      * data substance groups encoding HELMV2.0 polymer groups and extended
 *        annotations
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::string
rdkit_to_fasta(const ::RDKit::ROMol& mol);
} // namespace fasta
