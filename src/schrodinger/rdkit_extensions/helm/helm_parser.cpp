#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#include <algorithm>
#include <array>
#include <fmt/format.h>
#include <iterator>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"
#include "schrodinger/rdkit_extensions/helm/validation.h"

static void trim_string_view(std::string_view& value)
{
    while (!value.empty()) {
        if (std::isspace(value.front())) {
            value.remove_prefix(1);
        } else if (std::isspace(value.back())) {
            value.remove_suffix(1);
        } else {
            break;
        }
    }
};

[[nodiscard]] static std::array<std::string_view, 3>
parse_extended_annotations_and_version(const std::string_view input_helm)
{
    if (std::ranges::count(input_helm, '$') < 2) {
        return {input_helm, "", ""};
    }

    auto ext_annotations_end = input_helm.rfind('$');
    auto ext_annotations_start = input_helm.rfind('$', ext_annotations_end - 1);

    auto helm_prefix = input_helm.substr(0, ext_annotations_start + 1);
    auto extended_annotations =
        input_helm.substr(ext_annotations_start + 1,
                          ext_annotations_end - ext_annotations_start - 1);
    auto helm_version = ext_annotations_end == input_helm.size() - 1
                            ? ""
                            : input_helm.substr(ext_annotations_end + 1);

    return {helm_prefix, extended_annotations, helm_version};
}

namespace helm
{

HelmParser::HelmParser(const std::string_view input_helm) : m_input(input_helm)
{
    // strip surrounding whitespace
    trim_string_view(m_input);
}

helm_info HelmParser::parse()
{
    auto [input_helm_prefix, extended_annotations, helm_version] =
        parse_extended_annotations_and_version(m_input);

    std::istringstream stream(
        std::string{input_helm_prefix.data(), input_helm_prefix.size()});
    helm::TokenScanner scanner(&stream, m_input);
    helm::TokenParser parser(scanner, *this);

    // success return 0
    if (parser.parse() != 0) {
        std::stringstream ss;
        std::copy(m_errors.begin(), m_errors.end(),
                  std::ostream_iterator<std::string>(ss, "\n"));
        throw std::invalid_argument(ss.str());
    }

    addExtendedAnnotations(extended_annotations);
    addHelmVersion(helm_version);

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

void HelmParser::addExtendedAnnotations(
    const std::string_view extended_annotations)
{
    m_parsed_info.extended_annotations = extended_annotations;
}

void HelmParser::addHelmVersion(const std::string_view helm_version)
{
    m_parsed_info.helm_version = helm_version;
}

void HelmParser::addPolymerGroup(const std::string_view polymer_group_id,
                                 const std::string_view items,
                                 const bool is_polymer_union)
{
    m_parsed_info.polymer_groups.push_back(
        {polymer_group_id, items, is_polymer_union});
}

void HelmParser::addConnection(connection&& connection)
{
    m_parsed_info.connections.push_back(connection);
}

void HelmParser::addPolymer(const std::string_view polymer_id,
                            const std::string_view annotation)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.id = polymer_id;
    polymer.annotation = annotation;
}

void HelmParser::addMonomer(const std::string_view monomer_id,
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

void HelmParser::addMonomerWithId(const std::string_view monomer_id,
                                  const std::string_view annotation)
{
    addMonomer(monomer_id, false, false, false, annotation);
}

void HelmParser::addSmilesMonomer(const std::string_view smiles,
                                  const std::string_view annotation)
{
    addMonomer(smiles, true, false, false, annotation);
}

void HelmParser::addMonomerList(const std::string_view monomer_list,
                                const std::string_view annotation)
{
    addMonomer(monomer_list, false, false, true, annotation);
}

void HelmParser::markBranchMonomer(const size_t branch_group_size)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.monomers[polymer.monomers.size() - branch_group_size + 1]
        .is_branch = true;
}

void HelmParser::markLastNMonomersAsRepeated(
    const size_t repetition_size, const std::string_view num_repetitions,
    const std::string_view annotation)
{
    auto& polymer = m_parsed_info.polymers.back();
    polymer.repetitions.push_back({polymer.monomers.size() - repetition_size,
                                   repetition_size, num_repetitions,
                                   annotation});
}

void HelmParser::addResidueName(const std::string_view monomer_id)
{
    auto& polymers = m_parsed_info.polymers;
    if (polymers.empty() || !polymers.back().id.empty()) {
        polymers.push_back({});
    }
    auto& polymer = m_parsed_info.polymers.back();
    polymer.residue_names.insert(monomer_id);
}
void HelmParser::addWildcardOrUnknownResidue()
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
        "Malformed HELM string: check for mistakes around position {}:\n"
        "{}\n"
        "{}^\n"
        "{}",
        num_chars_processed, truncate_input(m_input, num_chars_processed - 1),
        std::string(num_dashes, '-'), err_msg));
}

[[nodiscard]] bool is_smiles_monomer(const std::string_view& token)
{
    // Assume it's SMILES if it has whitespace
    if (token.find_first_of(" \t\n\f\r\v") != std::string_view::npos) {
        return true;
    }

    auto monomer = token.substr(1, token.size() - 2);
    // Assume it's SMILES if it has these rgroup chars
    // NOTE: token is still surrounded by [] tokens
    if (monomer.find_first_of("*:[]=+") != std::string_view::npos) {
        return true;
    }

    // Parse the SMILES literally
    static const RDKit::v2::SmilesParse::SmilesParserParams smi_opts{
        .sanitize = false, .removeHs = false, .replacements = {}};

    [[maybe_unused]] schrodinger::rdkit_extensions::CaptureRDErrorLog rdkit_log;
    const std::string smiles{monomer};
    return RDKit::v2::SmilesParse::MolFromSmiles(smiles, smi_opts) != nullptr;
}
} // namespace helm
