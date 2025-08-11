#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <memory>
#include <string>

namespace RDKit
{
class RWMol;
}

namespace fasta
{
/*
 * Converts a fasta input to a monomeric mol. The fasta input can contain
 * multiple sequences, and each sequence can be annotated.
 *
 * @param peptide_fasta: one or more fasta sequences
 * @throws std::invalid_argument: If the input is malformed
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::RWMol>
peptide_fasta_to_rdkit(const std::string& peptide_fasta);

/*
 * Converts a fasta input to a monomeric mol. The fasta input can contain
 * multiple sequences, and each sequence can be annotated.
 *
 * @param rna_fasta: one or more fasta sequences
 * @throws std::invalid_argument: If the input is malformed
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::RWMol>
rna_fasta_to_rdkit(const std::string& rna_fasta);

/*
 * Converts a fasta input to a monomeric mol. The fasta input can contain
 * multiple sequences, and each sequence can be annotated.
 *
 * @param dna_fasta: one or more fasta sequences
 * @throws std::invalid_argument: If the input is malformed
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::RWMol>
dna_fasta_to_rdkit(const std::string& dna_fasta);
} // namespace fasta
