#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace fasta
{
[[nodiscard]] std::optional<std::string>
get_fasta_to_helm_amino_acid(char one_letter_monomer);

[[nodiscard]] std::optional<std::string>
get_fasta_to_helm_nucleotide(char one_letter_monomer,
                             std::string_view sugar_helm_monomer);

[[nodiscard]] std::optional<char>
get_helm_to_fasta_amino_acid(const std::string& helm_amino_acid);

[[nodiscard]] std::optional<char>
get_helm_to_fasta_nucleotide(const std::string& helm_nucleotide);
} // namespace fasta
