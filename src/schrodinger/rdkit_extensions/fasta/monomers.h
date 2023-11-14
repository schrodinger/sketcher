#pragma once

#include <optional>
#include <string_view>

namespace fasta
{
[[nodiscard]] std::optional<std::string_view>
get_fasta_to_helm_amino_acid(char one_letter_monomer);

[[nodiscard]] std::optional<std::string_view>
get_fasta_to_helm_nucleotide(char one_letter_monomer,
                             std::string_view sugar_helm_monomer);

} // namespace fasta
