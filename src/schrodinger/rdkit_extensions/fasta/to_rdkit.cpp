#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <regex>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/fasta/monomers.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

namespace
{
struct [[nodiscard]] fasta_sequence {
    std::optional<std::string> annotation;
    std::string sequence;
};
} // namespace

[[nodiscard]] static std::vector<fasta_sequence>
get_fasta_sequences(const std::string& generic_fasta);

template <class T> static std::unique_ptr<::RDKit::RWMol> get_monomer_mol
    [[nodiscard]] (const std::vector<fasta_sequence>& fasta_sequences,
                   T get_helm_monomer_function,
                   std::string_view polymer_prefix);

template <class T> static std::string get_polymer_helm
    [[nodiscard]] (const fasta_sequence& generic_fasta,
                   T get_helm_monomer_function, std::string_view polymer_prefix,
                   int polymer_number);

[[nodiscard]] static std::string
get_error_message(std::string_view sequence, unsigned int num_chars_processed,
                  std::string_view err_msg);

namespace fasta
{
[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
peptide_fasta_to_rdkit(const std::string& peptide_fasta)
{
    static constexpr std::string_view peptide_polymer_prefix{"PEPTIDE"};

    return get_monomer_mol(get_fasta_sequences(peptide_fasta),
                           get_fasta_to_helm_amino_acid,
                           peptide_polymer_prefix);
};

[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
rna_fasta_to_rdkit(const std::string& rna_fasta)
{
    static constexpr std::string_view rna_polymer_prefix{"RNA"};
    static constexpr std::string_view sugar_helm_monomer{"R"};

    return get_monomer_mol(
        get_fasta_sequences(rna_fasta),
        [](char one_letter_monomer) {
            return get_fasta_to_helm_nucleotide(one_letter_monomer,
                                                sugar_helm_monomer);
        },
        rna_polymer_prefix);
}

[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
dna_fasta_to_rdkit(const std::string& dna_fasta)
{
    static constexpr std::string_view dna_polymer_prefix{"RNA"};
    static constexpr std::string_view sugar_helm_monomer{"[dR]"};

    return get_monomer_mol(
        get_fasta_sequences(dna_fasta),
        [](char one_letter_monomer) {
            return get_fasta_to_helm_nucleotide(one_letter_monomer,
                                                sugar_helm_monomer);
        },
        dna_polymer_prefix);
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

template <class T> static std::string get_polymer_helm
    [[nodiscard]] (const fasta_sequence& generic_fasta,
                   T get_helm_monomer_function, std::string_view polymer_prefix,
                   int polymer_number)
{
    using namespace fmt::literals;

    std::vector<std::string> monomers;
    monomers.reserve(generic_fasta.sequence.size());
    std::transform(
        generic_fasta.sequence.begin(), generic_fasta.sequence.end(),
        std::back_inserter(monomers), [&, i = 0u](auto monomer) mutable {
            ++i;
            if (auto helm_monomer = get_helm_monomer_function(monomer);
                helm_monomer != std::nullopt) {
                return *helm_monomer;
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

template <class T> static std::unique_ptr<::RDKit::RWMol> get_monomer_mol
    [[nodiscard]] (const std::vector<fasta_sequence>& fasta_sequences,
                   T get_helm_monomer_function, std::string_view polymer_prefix)
{
    using namespace fmt::literals;

    std::vector<std::string> polymers;
    std::transform(
        fasta_sequences.begin(), fasta_sequences.end(),
        std::back_inserter(polymers), [&, i = 0](const auto& sequence) mutable {
            return get_polymer_helm(sequence, get_helm_monomer_function,
                                    polymer_prefix, ++i);
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

[[nodiscard]] static std::string
get_error_message(std::string_view sequence, unsigned int num_chars_processed,
                  std::string_view err_msg)
{
    // NOTE: If the input is very long, the pointer to the failed location
    // becomes less useful. We should truncate the length of the error message
    // to 101 chars.
    static constexpr int error_size{101};
    static constexpr unsigned int prefix_size{error_size / 2};
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
