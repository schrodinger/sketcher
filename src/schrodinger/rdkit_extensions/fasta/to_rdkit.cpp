#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fmt/format.h>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <regex>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

namespace
{
struct [[nodiscard]] fasta_sequence {
    std::optional<std::string> annotation;
    std::string sequence;
};
} // namespace

[[nodiscard]] static std::unordered_map<char, std::string>
get_fasta_to_helm_nucleotide_map(std::string_view sugar_helm_monomer);

[[nodiscard]] static std::vector<fasta_sequence>
get_fasta_sequences(const std::string& generic_fasta);

[[nodiscard]] static std::unique_ptr<::RDKit::RWMol> get_coarse_grain_rwmol(
    const std::vector<fasta_sequence>& fasta_sequences,
    const std::unordered_map<char, std::string>& fasta_to_helm_map,
    std::string_view polymer_prefix, std::string_view unknown_monomer_id);

[[nodiscard]] static std::string
get_polymer_helm(const fasta_sequence& generic_fasta,
                 const std::unordered_map<char, std::string>& fasta_to_helm_map,
                 std::string_view polymer_prefix, int polymer_number,
                 std::string_view unknown_monomer_id);

[[nodiscard]] static std::string get_error_message(std::string_view sequence,
                                                   int num_chars_processed,
                                                   std::string_view err_msg);

namespace fasta
{
[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
peptide_fasta_to_rdkit(const std::string& peptide_fasta)
{
    static const std::unordered_map<char, std::string>
        fasta_to_helm_amino_acid_map{
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
        };

    static constexpr std::string_view peptide_polymer_prefix{"PEPTIDE"};
    static constexpr std::string_view unknown_amino_acid{"X"};

    return get_coarse_grain_rwmol(get_fasta_sequences(peptide_fasta),
                                  fasta_to_helm_amino_acid_map,
                                  peptide_polymer_prefix, unknown_amino_acid);
};

[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
rna_fasta_to_rdkit(const std::string& rna_fasta)
{
    static constexpr std::string_view rna_polymer_prefix{"RNA"};
    static constexpr std::string_view unknown_nucleotide{"N"};
    static constexpr std::string_view sugar_helm_monomer{"R"};

    static const auto fasta_to_helm_rna_nucleotide_map =
        get_fasta_to_helm_nucleotide_map(sugar_helm_monomer);
    return get_coarse_grain_rwmol(get_fasta_sequences(rna_fasta),
                                  fasta_to_helm_rna_nucleotide_map,
                                  rna_polymer_prefix, unknown_nucleotide);
}

[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
dna_fasta_to_rdkit(const std::string& dna_fasta)
{
    static constexpr std::string_view dna_polymer_prefix{"RNA"};
    static constexpr std::string_view unknown_nucleotide{"N"};
    static constexpr std::string_view sugar_helm_monomer{"[dR]"};

    static const auto fasta_to_helm_dna_nucleotide_map =
        get_fasta_to_helm_nucleotide_map(sugar_helm_monomer);
    return get_coarse_grain_rwmol(get_fasta_sequences(dna_fasta),
                                  fasta_to_helm_dna_nucleotide_map,
                                  dna_polymer_prefix, unknown_nucleotide);
}
} // namespace fasta

// parses fasta sequences from the input text. some notes:
//  * recognizes '>' prefixed descriptions as sequence annotations
//  * sequences must immediately follow a sequence annotation, else they'll be
//    considered as a part of the previous sequence, if there is one
//  * recognizes sequences/subsequences separated by whitespace as a single
//    sequence. i.e.
//                  ABC
//                  DEF
//    and
//                  A B
//                      C D
//                  EF
//
//    will be recognized as ABCDEF
//  * sequences can contain gaps and spaces -  these will not be functionally
//    significant.
//  * translation stops '*' are not supported
[[nodiscard]] static std::vector<fasta_sequence>
get_fasta_sequences(const std::string& generic_fasta)
{
    std::vector<fasta_sequence> fasta_sequences;
    std::optional<fasta_sequence> current_sequence;

    auto save_current_sequence = [&]() {
        if (current_sequence && current_sequence->annotation) {
            fasta_sequences.push_back(std::move(*current_sequence));
            current_sequence = std::nullopt;
        }
    };

    std::vector<std::string> lines;
    boost::algorithm::split(lines, generic_fasta + "\n>",
                            boost::is_any_of("\n"));
    for (const auto& line : lines) {
        if (line.empty()) {
            continue; // skip line
        } else if (line[0] == '>') {
            save_current_sequence();
            current_sequence = {line.substr(1), ""};
        } else if (current_sequence) {
            current_sequence->sequence.append(std::move(line));
        } else {
            current_sequence = {std::nullopt, std::move(line)};
        }
    }

    static const std::regex chars_to_delete{R"([- \n])"};
    std::for_each(
        fasta_sequences.begin(), fasta_sequences.end(),
        [&](auto& fasta_sequence) {
            auto& [annotation, monomers] = fasta_sequence;
            monomers = std::regex_replace(monomers, chars_to_delete, "");

            if (monomers.empty() && annotation) {
                throw std::invalid_argument(fmt::format(
                    "There's no sequence associated with annotation, '{}'",
                    *annotation));
            }
        });
    return fasta_sequences;
}

[[nodiscard]] static std::string
get_polymer_helm(const fasta_sequence& generic_fasta,
                 const std::unordered_map<char, std::string>& fasta_to_helm_map,
                 std::string_view polymer_prefix, int polymer_number,
                 std::string_view unknown_monomer_id)
{
    using namespace fmt::literals;

    std::vector<std::string_view> monomers;
    monomers.reserve(generic_fasta.sequence.size());
    std::transform(generic_fasta.sequence.begin(), generic_fasta.sequence.end(),
                   std::back_inserter(monomers),
                   [&, i = 0](auto monomer) mutable -> std::string_view {
                       ++i;
                       if (fasta_to_helm_map.count(monomer)) {
                           return fasta_to_helm_map.at(monomer);
                       }

                       // char vs string_view
                       if (monomer == unknown_monomer_id[0]) {
                           return unknown_monomer_id;
                       }

                       throw std::invalid_argument(get_error_message(
                           generic_fasta.sequence, i, "Unsupported monomer"));
                   });

    static constexpr std::string_view monomer_separator{"."};
    auto monomer_sequence = fmt::join(monomers, monomer_separator);

    if (generic_fasta.annotation && !generic_fasta.annotation->empty()) {
        static constexpr std::string_view helm_polymer_template{
            "{polymer_prefix}{polymer_number}{{{monomer_sequence}}}\"{"
            "annotation}\""};
        return fmt::format(helm_polymer_template,
                           "polymer_prefix"_a = polymer_prefix,
                           "polymer_number"_a = polymer_number,
                           "monomer_sequence"_a = monomer_sequence,
                           "annotation"_a = *generic_fasta.annotation);
    } else {
        static constexpr std::string_view helm_polymer_template{
            "{polymer_prefix}{polymer_number}{{{monomer_sequence}}}"};
        return fmt::format(helm_polymer_template,
                           "polymer_prefix"_a = polymer_prefix,
                           "polymer_number"_a = polymer_number,
                           "monomer_sequence"_a = monomer_sequence);
    }
}

[[nodiscard]] static std::unique_ptr<::RDKit::RWMol> get_coarse_grain_rwmol(
    const std::vector<fasta_sequence>& fasta_sequences,
    const std::unordered_map<char, std::string>& fasta_to_helm_map,
    std::string_view polymer_prefix, std::string_view unknown_monomer_id)
{
    using namespace fmt::literals;

    std::vector<std::string> polymers;
    std::transform(
        fasta_sequences.begin(), fasta_sequences.end(),
        std::back_inserter(polymers), [&, i = 0](const auto& sequence) mutable {
            return get_polymer_helm(sequence, fasta_to_helm_map, polymer_prefix,
                                    ++i, unknown_monomer_id);
        });

    if (polymers.empty()) {
        throw std::invalid_argument(
            "Couldn't retrieve any FASTA sequences from input");
    }

    try {
        static constexpr std::string_view polymer_separator{"|"};

        return helm::helm_to_rdkit(
            fmt::format("{polymers}$$$$V2.0",
                        "polymers"_a = fmt::join(polymers, polymer_separator)));
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Couldn't convert FASTA sequence to mol.");
    }
}

[[nodiscard]] static std::string get_error_message(std::string_view sequence,
                                                   int num_chars_processed,
                                                   std::string_view err_msg)
{
    // NOTE: If the input is very long, the pointer to the failed location
    // becomes less useful. We should truncate the length of the error message
    // to 101 chars.
    static constexpr int error_size{101};
    static constexpr int prefix_size{error_size / 2};
    static auto truncate_input = [](const auto& input, const unsigned int pos) {
        if ((pos >= prefix_size) && (pos + prefix_size) < input.size()) {
            return input.substr(pos - prefix_size, error_size);
        } else if (pos >= prefix_size) {
            return input.substr(pos - prefix_size);
        } else {
            return input.substr(
                0, std::min(input.size(), static_cast<size_t>(error_size)));
        }
    };

    size_t num_dashes =
        (num_chars_processed >= prefix_size ? prefix_size
                                            : num_chars_processed - 1);
    return fmt::format(
        "Malformed FASTA string: check for mistakes around position {}:\n"
        "{}\n"
        "{}^\n"
        "{}",
        num_chars_processed, truncate_input(sequence, num_chars_processed - 1),
        std::string(num_dashes, '-'), err_msg);
}

[[nodiscard]] static std::unordered_map<char, std::string>
get_fasta_to_helm_nucleotide_map(std::string_view sugar_helm_monomer)
{
    // NOTE: the helm entries have to separate the phosphate, sugar and
    // nucleotide base subunits
    static const std::unordered_map<char, std::string_view>
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
