#include "schrodinger/rdkit_extensions/fasta/monomers.h"

#include <algorithm>
#include <fmt/format.h>
#include <initializer_list>

namespace fasta
{

[[nodiscard]] static std::unordered_map<char, std::string_view>
get_fasta_to_helm_amino_acid_map()
{
    return {{
        {'A', "A"}, // alanine
        {'C', "C"}, // cystine
        {'D', "D"}, // aspartate
        {'E', "E"}, // glutamate
        {'F', "F"}, // phenylalanine
        {'G', "G"}, // glycine
        {'H', "H"}, // histidine
        {'I', "I"}, // isoleucine
        {'K', "K"}, // lysine
        {'L', "L"}, // leucine
        {'M', "M"}, // methionine
        {'N', "N"}, // asparagine
        {'P', "P"}, // proline
        {'Q', "Q"}, // glutamine
        {'R', "R"}, // arginine
        {'S', "S"}, // serine
        {'T', "T"}, // threonine
        {'V', "V"}, // valine
        {'W', "W"}, // tryptophan
        {'Y', "Y"}, // tyrosine
        // Non-standard
        {'O', "O"}, // pyrrolysine
        {'U', "U"}, // selenocysteine
        // Ambiguous
        {'B', "(D+N)"}, // aspartate or asparagine
        {'J', "(I+L)"}, // isoleucine or leucine
        {'Z', "(E+Q)"}, // glutamate or glutamine
        // Unknown monomer
        {'X', "X"}, // any, meaning one single unknown amino acid
    }};
};

[[nodiscard]] static std::unordered_map<char, std::string>
get_fasta_to_helm_nucleotide_map(std::string_view sugar_helm_monomer)
{
    // NOTE: the helm entries have to separate the phosphate, sugar and
    // nucleotide base subunits
    static constexpr std::initializer_list<std::pair<char, std::string_view>>
        fasta_to_helm_nucleotide_template_map{
            {'A', "{}(A)P"}, // adenosine
            {'C', "{}(C)P"}, // cytidine
            {'G', "{}(G)P"}, // guanine
            {'U', "{}(U)P"}, // uridine
            {'T', "{}(T)P"}, // thymidine
            // Ambiguous
            {'K', "{}((G+T))P"},   // keto
            {'S', "{}((G+C))P"},   // strong
            {'Y', "{}((T+C))P"},   // pyrimidine
            {'M', "{}((A+C))P"},   // amino
            {'W', "{}((A+T))P"},   // weak
            {'R', "{}((G+A))P"},   // purine
            {'B', "{}((G+T+C))P"}, // G/T/C
            {'D', "{}((G+A+T))P"}, // G/A/T
            {'H', "{}((A+C+T))P"}, // A/C/T
            {'V', "{}((G+C+A))P"}, // G/C/A
            // Unknown monomer
            {'N', "{}(N)P"}, // any, meaning one single unknown base
        };

    std::unordered_map<char, std::string> nucleotide_map;
    std::transform(
        fasta_to_helm_nucleotide_template_map.begin(),
        fasta_to_helm_nucleotide_template_map.end(),
        std::inserter(nucleotide_map, nucleotide_map.end()),
        [&](const auto& fasta_to_helm_pair) {
            auto& [fasta_nucleotide, helm_template] = fasta_to_helm_pair;
            return std::make_pair(
                fasta_nucleotide,
                fmt::format(fmt::runtime(helm_template), sugar_helm_monomer));
        });
    return nucleotide_map;
};

[[nodiscard]] std::optional<std::string_view>
get_fasta_to_helm_amino_acid(char one_letter_monomer)
{
    static const auto fasta_to_helm_map = get_fasta_to_helm_amino_acid_map();

    auto helm_monomer = fasta_to_helm_map.find(one_letter_monomer);
    return helm_monomer == fasta_to_helm_map.end()
               ? std::nullopt
               : std::optional<std::string_view>(helm_monomer->second);
};

[[nodiscard]] std::optional<std::string_view>
get_fasta_to_helm_nucleotide(char one_letter_monomer,
                             std::string_view sugar_helm_monomer)
{
    static const auto rna_map = get_fasta_to_helm_nucleotide_map("R");
    static const auto dna_map = get_fasta_to_helm_nucleotide_map("[dR]");

    if (sugar_helm_monomer != "R" && sugar_helm_monomer != "[dR]") {
        throw std::invalid_argument(
            "FASTA conversions with sugar monomers that aren't 'R' or '[dR]' "
            "is currently unsupported.");
    }

    auto& fasta_to_helm_map = sugar_helm_monomer == "R" ? rna_map : dna_map;
    auto helm_monomer = fasta_to_helm_map.find(one_letter_monomer);
    return helm_monomer == fasta_to_helm_map.end()
               ? std::nullopt
               : std::optional<std::string_view>(helm_monomer->second);
}
} // namespace fasta
