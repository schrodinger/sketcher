#include "schrodinger/rdkit_extensions/helm/validation.h"

#include <algorithm>
#include <boost/json.hpp>
#include <charconv> // from_chars
#include <cctype>
#include <fmt/format.h>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

namespace helm
{
namespace
{
[[nodiscard]] bool validate_smiles(const std::string_view& monomer_smiles,
                                   const size_t residue_number,
                                   const size_t num_monomers,
                                   const bool is_branch_monomer)
{
    // NOTE: Types of Rgroups in smiles, [OH:1], [*:1]
    static std::regex modern_rgroup_pattern(R"(\[\w*\*?\:([1-9][0-9]*)\])");
    static std::regex cxxsmiles_rgroup_pattern(R"(_?R([1-9][0-9]*))");

    // Make sure we inline smiles have the attachment points needed to
    // for the right bonds
    auto check_smiles = [&](const std::regex& rgroup_pattern) {
        std::string smiles{monomer_smiles.begin(), monomer_smiles.end()};
        std::sregex_token_iterator matches(smiles.begin(), smiles.end(),
                                           rgroup_pattern, 1);

        std::unordered_set<unsigned int> rgroup_ids;
        std::transform(matches, std::sregex_token_iterator(),
                       std::inserter(rgroup_ids, rgroup_ids.end()),
                       [](const auto& id) { return std::stoi(id.str()); });

        // If this is the only residue in the polymer, we don't need to look for
        // rgroups now
        if (residue_number == 1 && num_monomers == 1) {
            return true;
        }
        // Anything else should have rgroups
        if (!rgroup_ids.size()) {
            return false;
        }
        // We need R1 to form the bond to previous atom. i.e. backbone bond
        // would be R2-R1 and branch bond would be R3-R1
        if (residue_number > 1 && !rgroup_ids.count(1)) {
            return false;
        }
        // We need either or both of R2 and R3 to form the next bond if this
        // isn't a branch monomer
        if (residue_number < num_monomers && !is_branch_monomer &&
            !(rgroup_ids.count(2) || rgroup_ids.count(3))) {
            return false;
        }
        return true;
    };

    return check_smiles(modern_rgroup_pattern) ||
           check_smiles(cxxsmiles_rgroup_pattern);
}

void validate_polymers(const std::vector<polymer>& polymers,
                       std::unordered_set<std::string_view>& polymer_ids,
                       std::vector<std::string>& errors,
                       const std::string_view& input)
{
    for (const auto& polymer : polymers) {
        if (!polymer_ids.insert(polymer.id).second) {
            auto msg = fmt::format("Duplicate polymer id: '{}'", polymer.id);
            errors.push_back(construct_error_msg(msg, polymer.id, input));
        }

        const auto num_monomers = polymer.monomers.size();
        auto i = 0;
        for (const auto& monomer : polymer.monomers) {
            ++i;
            if (monomer.is_smiles &&
                !validate_smiles(monomer.id, i, num_monomers,
                                 monomer.is_branch)) {
                auto msg = "The inline smiles has missing rgroups.";
                errors.push_back(construct_error_msg(msg, monomer.id, input));
            } else if (monomer.id == "_") {
                auto msg = "Missing monomers are currently unsupported.";
                errors.push_back(construct_error_msg(msg, monomer.id, input));
            } else if (monomer.is_list &&
                       monomer.id.find("_") != std::string_view::npos) {
                auto token = monomer.id.substr(monomer.id.find("_"));
                auto msg = "Missing monomers are currently unsupported.";
                errors.push_back(construct_error_msg(msg, token, input));
            }
        }
    }
}

void validate_connection_info(
    const std::string_view& polymer_id, const std::string_view& monomer_id,
    const std::string_view& rgroup,
    const std::unordered_set<std::string_view>& polymer_ids,
    std::vector<std::string>& errors, const std::string_view& input)
{
    // polymer ids must exist
    if (polymer_id[0] != 'G' && !polymer_ids.count(polymer_id)) {
        auto msg =
            fmt::format("Unknown polymer '{}' in connection", polymer_id);
        errors.push_back(construct_error_msg(msg, polymer_id, input));
    }

    if (polymer_id.starts_with("BLOB")) {
        if (monomer_id != "*" && monomer_id != "?") {
            auto msg = "Connections to/from BLOB polymers  can't have "
                       "defined connection residues";
            errors.push_back(construct_error_msg(msg, monomer_id, input));
        }

        if (rgroup != "*" && rgroup != "?") {
            auto msg = "Connections to/from BLOB polymers  can't have "
                       "defined connection rgroups";
            errors.push_back(construct_error_msg(msg, rgroup, input));
        }
    } else {
        if (monomer_id.find_first_of("*?") != std::string_view::npos) {
            auto msg = "Ambiguous bonds to non-blob polymers are currently "
                       "unsupported";
            errors.push_back(construct_error_msg(msg, monomer_id, input));
        }

        if (rgroup.find_first_of("*?") != std::string_view::npos) {
            auto msg = "Ambiguous bonds to non-blob polymers are currently "
                       "unsupported";
            errors.push_back(construct_error_msg(msg, rgroup, input));
        }
    }

    // Unsupported features
    if (polymer_id[0] == 'G') {
        auto msg = "Connections involving polymer groups are currently "
                   "unsupported.";
        errors.push_back(construct_error_msg(msg, polymer_id, input));
    }

    if (polymer_id.find_first_of(",+") != std::string_view::npos) {
        auto msg = "Multiple monomers in connections are currently unsupported";
        errors.push_back(construct_error_msg(msg, polymer_id, input));
    }
}

//
// Rules for a valid connection:
//      * BLOBS can't have defined connection residues or rgroups
//      * Unknown residues can't have defined rgroups
//
// Unsupported connection types:
//      * Connections with polymer groups
//      * Connections with non-blob ambiguous bonds
//      * Connections with multiple to and from residues
void validate_connection(
    const connection& connection,
    const std::unordered_set<std::string_view>& polymer_ids,
    std::vector<std::string>& errors, const std::string_view& input)
{
    // check hydrogen pairing must be the same
    if ((connection.from_rgroup == "pair") ^ (connection.to_rgroup == "pair")) {
        auto msg = "Connections for hydrogen pairings must have 'pair' as "
                   "the attachment point for both monomers.";
        errors.push_back(construct_error_msg(msg, connection.from_id, input));
    }

    validate_connection_info(connection.from_id, connection.from_res,
                             connection.from_rgroup, polymer_ids, errors,
                             input);

    validate_connection_info(connection.to_id, connection.to_res,
                             connection.to_rgroup, polymer_ids, errors, input);
}

[[nodiscard]] std::vector<std::string_view>
split_string_view(const std::string_view& source,
                  const std::string_view delimiters)
{
    std::vector<std::string_view> tokens;
    size_t start = 0;
    for (size_t idx = 1; idx < source.size(); ++idx) {
        if (delimiters.find(source[idx]) != std::string_view::npos) {
            tokens.push_back(source.substr(start, idx - start));
            start = idx + 1;
        }
    }
    tokens.push_back(source.substr(start));
    return tokens;
}

void validate_connection_residues(const std::string_view& residue,
                                  const std::string_view& rgroup,
                                  const polymer& polymer,
                                  std::vector<std::string>& errors,
                                  const std::string_view& input)
{
    for (const auto& res : split_string_view(residue, "+,")) {
        size_t residue_number = 0;
        auto status = std::from_chars(res.data(), res.data() + res.size(),
                                      residue_number);
        if (status.ec == std::errc()) { // residue number
            if (residue_number > polymer.monomers.size()) {
                auto msg = fmt::format(
                    "Residue number '{}' falls outside of polymer chain.",
                    residue_number);
                errors.push_back(construct_error_msg(msg, res, input));
            } else if (polymer.wildcard_and_unknown_residues.count(
                           residue_number) &&
                       rgroup[0] == 'R') {
                auto msg = fmt::format("Residue number '{}' points to a "
                                       "wildcard/ unknown residue, so "
                                       "it can't have a specified rgroup.",
                                       residue_number);
                errors.push_back(construct_error_msg(msg, res, input));
            }
        } else if (!polymer.residue_names.count(residue)) {
            auto msg =
                fmt::format("Residue '{}' can't be found in polymer.", res);
            errors.push_back(construct_error_msg(msg, res, input));
        }
    }
}

void validate_non_blob_connection_residues(const std::string_view& polymer_id,
                                           const std::string_view& residue,
                                           const std::string_view& rgroup,
                                           const std::vector<polymer>& polymers,
                                           std::vector<std::string>& errors,
                                           const std::string_view& input)
{
    // non-blob polymers should have reachable residues
    auto polymer = std::ranges::find_if(polymers, [&](const auto& polymer) {
        return polymer.id == polymer_id;
    });

    if (polymer == polymers.end()) {
        auto msg =
            fmt::format("Unknown polymer '{}' in connection", polymer_id);
        errors.push_back(construct_error_msg(msg, polymer_id, input));
    } else {
        validate_connection_residues(residue, rgroup, *polymer, errors, input);
    }
}

// Helper api to validate all parsed connections
void validate_connections(
    const std::vector<connection>& connections,
    const std::vector<polymer>& polymers,
    const std::unordered_set<std::string_view>& polymer_ids,
    std::vector<std::string>& errors, const std::string_view& input)
{
    for (const auto& connection : connections) {
        auto& [from_id, from_res, from_rgroup, to_id, to_res, to_rgroup,
               unused] = connection;

        validate_connection(connection, polymer_ids, errors, input);

        // non-blob polymers should have reachable residues
        if (!connection.from_id.starts_with("BLOB")) {
            validate_non_blob_connection_residues(
                from_id, from_res, from_rgroup, polymers, errors, input);
        }

        // non-blob polymers should have reachable residues
        if (!connection.to_id.starts_with("BLOB")) {
            validate_non_blob_connection_residues(to_id, to_res, to_rgroup,
                                                  polymers, errors, input);
        }
    }
}

void validate_polymer_groups(
    const std::vector<polymer_group>& polymer_groups,
    const std::unordered_set<std::string_view>& polymer_ids,
    std::vector<std::string>& errors, const std::string_view& input)
{
    std::unordered_set<std::string_view> group_ids;
    // check duplicate group ids
    for (const auto& polymer_group : polymer_groups) {
        if (!group_ids.insert(polymer_group.id).second) {
            auto msg = fmt::format("Duplicate polymer group id: '{}'",
                                   polymer_group.id);
            errors.push_back(construct_error_msg(msg, polymer_group.id, input));
        }
    }

    // Make sure all items in the list can be located
    for (const auto& polymer_group : polymer_groups) {
        for (const auto& item : split_string_view(polymer_group.items, "+,")) {
            std::string_view prefix = item;
            // strip ratio off
            if (prefix.find(":") != std::string_view::npos) {
                prefix.remove_suffix(prefix.size() - prefix.find(":"));
            }

            if (!(polymer_ids.count(prefix) || group_ids.count(prefix))) {
                auto msg =
                    fmt::format("Unknown polymer group item: '{}'", item);
                errors.push_back(construct_error_msg(msg, item, input));
            }
        }
    }
}

void validate_extended_annotations(std::string_view extended_annotations,
                                   std::vector<std::string>& errors,
                                   const std::string_view& input)
{
    if (!extended_annotations.empty()) {
        boost::json::stream_parser p;
        boost::system::error_code ec;
        p.write(extended_annotations.data(), extended_annotations.size(), ec);
        if (ec) {
            auto msg = "Invalid extended annotations";
            errors.push_back(
                construct_error_msg(msg, extended_annotations, input));
        }
    }
}

void validate_helm_version(const std::string_view& helm_version,
                           const std::string_view& extended_annotations,
                           std::vector<std::string>& errors,
                           const std::string_view& input)
{
    if (!helm_version.empty() && helm_version != "V2.0") {
        auto msg = "Only HELM and HELMV2.0 versions are currently supported";
        errors.push_back(construct_error_msg(msg, helm_version, input));
    }

    if (helm_version.empty() && !extended_annotations.empty()) {
        auto msg = "HELM annotations are currently unsupported";
        errors.push_back(construct_error_msg(msg, extended_annotations, input));
    }
}
} // namespace

bool validate_parsed_info(const helm_info& parsed_info,
                          std::vector<std::string>& errors,
                          const std::string_view& input)
{
    auto& [polymers, connections, polymer_groups, extended_annotations,
           helm_version] = parsed_info;

    std::unordered_set<std::string_view> polymer_ids;
    validate_polymers(polymers, polymer_ids, errors, input);
    validate_connections(connections, polymers, polymer_ids, errors, input);
    validate_polymer_groups(polymer_groups, polymer_ids, errors, input);
    validate_extended_annotations(extended_annotations, errors, input);
    validate_helm_version(helm_version, extended_annotations, errors, input);

    return errors.empty();
}

} // namespace helm
