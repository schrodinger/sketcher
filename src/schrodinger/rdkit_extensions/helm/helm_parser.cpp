#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fmt/format.h>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"
#include "schrodinger/rdkit_extensions/helm/validation.h"

namespace helm
{

helm_info HelmParser::parse()
{
    std::istringstream stream(std::string{m_input});
    helm::TokenScanner scanner(&stream, m_input);
    helm::TokenParser parser(scanner, *this);

    // success return 0
    if (parser.parse() != 0) {
        std::stringstream ss;
        std::copy(m_errors.begin(), m_errors.end(),
                  std::ostream_iterator<std::string>(ss, "\n"));
        throw std::invalid_argument(ss.str());
    }

    validate_parsed_info(m_parsed_info, *this);
    if (hasErrors()) {
        std::stringstream ss;
        ss << "Parsing failed because of the following: \n";
        std::copy(m_errors.begin(), m_errors.end(),
                  std::ostream_iterator<std::string>(ss, "\n"));
        throw std::invalid_argument(ss.str());
    }

    return std::move(m_parsed_info);
}

void HelmParser::add_extended_annotation(
    const std::string_view extended_annotation)
{
    m_parsed_info.extended_annotation = extended_annotation;
}

void HelmParser::add_polymer_group(const std::string_view polymer_group_id,
                                   const std::string_view items,
                                   const bool is_polymer_union)
{
    m_parsed_info.polymer_groups.push_back(
        {polymer_group_id, items, is_polymer_union});
}

void HelmParser::add_connection(connection&& connection)
{
    m_parsed_info.connections.push_back(connection);
}

void HelmParser::add_polymer(const std::string_view polymer_id,
                             const std::string_view annotation)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.id = polymer_id;
    polymer.annotation = annotation;
}

void HelmParser::add_monomer(const std::string_view monomer_id,
                             const bool is_smiles, const bool is_branch,
                             const bool is_list,
                             const std::string_view annotation)
{
    auto& polymers = m_parsed_info.polymers;
    if (polymers.empty() || !polymers.back().id.empty()) {
        polymers.push_back({});
    }
    auto& polymer = polymers.back();
    polymer.monomers.push_back(
        {monomer_id, is_smiles, is_branch, is_list, annotation});
}

void HelmParser::add_monomer_with_id(const std::string_view monomer_id,
                                     const std::string_view annotation)
{
    add_monomer(monomer_id, false, false, false, annotation);
}

void HelmParser::add_smiles_monomer(const std::string_view smiles,
                                    const std::string_view annotation)
{
    add_monomer(smiles, true, false, false, annotation);
}

void HelmParser::add_monomer_list(const std::string_view monomer_list,
                                  const std::string_view annotation)
{
    add_monomer(monomer_list, false, false, true, annotation);
}

void HelmParser::mark_branch_monomer(const size_t branch_group_size)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.monomers[polymer.monomers.size() - branch_group_size + 1]
        .is_branch = true;
}

void HelmParser::mark_last_n_monomers_as_repeated(
    const size_t repetition_size, const std::string_view num_repetitions,
    const std::string_view annotation)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.repetitions.push_back({polymer.monomers.size() - repetition_size,
                                   repetition_size, num_repetitions,
                                   annotation});
}

void HelmParser::add_residue_name(const std::string_view monomer_id)
{
    auto& polymers = m_parsed_info.polymers;
    if (polymers.empty() || !polymers.back().id.empty()) {
        polymers.push_back({});
    }
    auto& polymer = m_parsed_info.polymers.back();
    polymer.residue_names.insert(monomer_id);
}
void HelmParser::add_wildcard_or_unknown_residue()
{
    auto& polymers = m_parsed_info.polymers;
    if (polymers.empty() || !polymers.back().id.empty()) {
        polymers.push_back({});
    }
    auto& polymer = m_parsed_info.polymers.back();
    polymer.wildcard_and_unknown_residues.insert(polymer.monomers.size() + 1);
}

bool HelmParser::hasErrors()
{
    return !m_errors.empty();
}

void HelmParser::saveError(const std::string_view& failed_token,
                           const std::string& err_msg)
{
    const auto num_chars_processed =
        std::distance(m_input.data(), failed_token.data()) + 1;
    saveError(num_chars_processed, err_msg);
}

void HelmParser::saveError(const unsigned int num_chars_processed,
                           const std::string& err_msg)
{
    // NOTE: If the input is very long, the pointer to the failed location
    // becomes less useful. We should truncate the length of the error message
    // to 101 chars.
    static constexpr unsigned int error_size{101};
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
    m_errors.push_back(fmt::format(
        "Parsing failed around position {}:\n"
        "{}\n"
        "{}^\n"
        "{}",
        num_chars_processed, truncate_input(m_input, num_chars_processed - 1),
        std::string(num_dashes, '-'), err_msg));
}

} // namespace helm
