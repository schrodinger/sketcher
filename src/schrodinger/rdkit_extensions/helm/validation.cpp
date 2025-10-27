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

void validate_polymer(const polymer& polymer,
                      std::unordered_set<std::string_view>& polymer_ids,
                      HelmParser& parser)
{
    if (!polymer_ids.insert(polymer.id).second) {
        parser.saveError(polymer.id,
                         fmt::format("Duplicate polymer id: '{}'", polymer.id));
    }

    const auto num_monomers = polymer.monomers.size();
    auto i = 0;
    for (const auto& monomer : polymer.monomers) {
        ++i;
        if (monomer.is_smiles &&
            !validate_smiles(monomer.id, i, num_monomers, monomer.is_branch)) {
            parser.saveError(monomer.id,
                             "The inline smiles has missing rgroups.");
        } else if (monomer.id == "_") {
            parser.saveError(monomer.id,
                             "Missing monomers are currently unsupported.");
        } else if (monomer.is_list &&
                   monomer.id.find("_") != std::string_view::npos) {

            parser.saveError(monomer.id.substr(monomer.id.find("_")),
                             "Missing monomers are currently unsupported.");
        }
    }
}

/*
 * Rules for a valid connection:
 *      * BLOBS can't have defined connection residues or rgroups
 *      * Unknown residues can't have defined rgroups
 *
 * Unsupported connection types:
 *      * Connections with polymer groups
 *      * Connections with non-blob ambiguous bonds
 *      * Connections with multiple to and from residues
 */
void validate_connection(
    const connection& connection,
    const std::unordered_set<std::string_view>& polymer_ids, HelmParser& parser)
{
    static auto validate_blob_connection_info = [](const auto& monomer_id,
                                                   const auto& rgroup,
                                                   auto& parser) {
        if (monomer_id != "*" && monomer_id != "?") {
            parser.saveError(monomer_id,
                             "Connections to/from BLOB polymers  can't have "
                             "defined connection residues");
        }
        if (rgroup != "*" && rgroup != "?") {
            parser.saveError(rgroup,
                             "Connections to/from BLOB polymers  can't have "
                             "defined connection rgroups");
        }
    };

    static auto validate_hydrogen_pairings = [](const auto& connection,
                                                auto& parser) {
        if (connection.from_rgroup == "pair" &&
            connection.to_rgroup == "pair") {
            return;
        }

        if (connection.from_rgroup == "pair" ||
            connection.to_rgroup == "pair") {
            parser.saveError(
                connection.from_id,
                "Connections for hydrogen pairings must have 'pair' as the "
                "attachment point for both monomers.");
        }
    };

    static auto validate_connection_info = [](const auto& polymer_ids,
                                              const auto& polymer_id,
                                              const auto& monomer_id,
                                              const auto& rgroup,
                                              auto& parser) {
        if (polymer_id[0] != 'G' && !polymer_ids.count(polymer_id)) {
            parser.saveError(
                polymer_id,
                fmt::format("Unknown polymer '{}' in connection", polymer_id));
        } else if (polymer_id[0] == 'B') {
            validate_blob_connection_info(monomer_id, rgroup, parser);
        } else if (polymer_id[0] == 'G') {
            parser.saveError(polymer_id,
                             "Connections involving polymer groups are "
                             "currently unsupported.");
        } else if (monomer_id.find_first_of("*?") != std::string_view::npos ||
                   rgroup.find_first_of("*?") != std::string_view::npos) {
            parser.saveError(rgroup, "Ambiguous bonds to non-blob polymers are "
                                     "currently unsupported");
        } else if (monomer_id.find_first_of(",+") != std::string_view::npos) {
            parser.saveError(polymer_id, "Multiple monomers in connections are "
                                         "currently unsupported");
        }
    };

    validate_hydrogen_pairings(connection, parser);

    validate_connection_info(polymer_ids, connection.from_id,
                             connection.from_res, connection.from_rgroup,
                             parser);
    validate_connection_info(polymer_ids, connection.to_id, connection.to_res,
                             connection.to_rgroup, parser);
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

void validate_connection_residues(const connection& connection,
                                  const std::vector<polymer>& polymers,
                                  HelmParser& parser)
{
    static auto validate_residue_number = [](const auto& residue,
                                             const auto& polymer,
                                             const auto& rgroup, auto& parser) {
        size_t residue_number = 0;
        std::from_chars(residue.data(), residue.data() + residue.size(),
                        residue_number);
        if (residue_number > polymer.monomers.size()) {
            parser.saveError(
                residue,
                fmt::format(
                    "Residue number '{}' falls outside of polymer chain.",
                    residue_number));
        } else if (polymer.wildcard_and_unknown_residues.count(
                       residue_number) &&
                   rgroup[0] == 'R') {
            parser.saveError(residue,
                             fmt::format("Residue number '{}' points to a "
                                         "wildcard/ unknown residue, so "
                                         "it can't have a specified rgroup.",
                                         residue_number));
        }
    };

    static auto validate_connection_residue =
        [](const auto& polymer, const auto& residues, const auto& rgroup,
           auto& parser) {
            for (const auto& residue : split_string_view(residues, "+,")) {
                const auto is_digit = std::all_of(
                    residue.begin(), residue.end(),
                    [](const unsigned char c) { return std::isdigit(c); });
                if (is_digit) {
                    validate_residue_number(residue, polymer, rgroup, parser);
                } else if (!polymer.residue_names.count(residue)) {
                    parser.saveError(
                        residue,
                        fmt::format("Residue '{}' can't be found in polymer.",
                                    residue));
                }
            }
        };

    static auto find_polymer = [](const auto& polymer_id,
                                  const auto& polymers) {
        return *std::find_if(
            polymers.begin(), polymers.end(),
            [&](const auto& polymer) { return polymer.id == polymer_id; });
    };

    if (connection.from_id[0] != 'B') {
        validate_connection_residue(find_polymer(connection.from_id, polymers),
                                    connection.from_res, connection.from_rgroup,
                                    parser);
    }
    if (connection.to_id[0] != 'B') {
        validate_connection_residue(find_polymer(connection.to_id, polymers),
                                    connection.to_res, connection.to_rgroup,
                                    parser);
    }
}
// check duplicate group ids
void validate_polymer_group_ids(const polymer_group& polymer_group,
                                std::unordered_set<std::string_view>& group_ids,
                                HelmParser& parser)
{
    if (!group_ids.insert(polymer_group.id).second) {
        parser.saveError(
            polymer_group.id,
            fmt::format("Duplicate polymer group id: '{}'", polymer_group.id));
    }
}

// Make sure all items in the list can be located
void validate_polymer_group_items(
    const polymer_group& polymer_group,
    const std::unordered_set<std::string_view>& polymer_ids,
    const std::unordered_set<std::string_view>& group_ids, HelmParser& parser)
{
    for (const auto& item : split_string_view(polymer_group.items, "+,")) {
        std::string_view prefix = item;
        // strip ratio off
        if (prefix.find(":") != std::string_view::npos) {
            prefix.remove_suffix(prefix.size() - prefix.find(":"));
        }
        if (!(polymer_ids.count(prefix) || group_ids.count(prefix))) {
            parser.saveError(
                prefix, fmt::format("Unknown polymer group item: '{}'", item));
        }
    }
}

void validate_extended_annotations(std::string_view extended_annotations,
                                   HelmParser& parser)
{
    if (extended_annotations.empty()) {
        return;
    }

    boost::json::stream_parser p;
    boost::system::error_code ec;
    p.write(extended_annotations.data(), extended_annotations.size(), ec);
    if (ec) {
        parser.saveError(extended_annotations, "Invalid extended annotations");
    }
}

void validate_helm_version(const std::string_view& helm_version,
                           const std::string_view& extended_annotations,
                           HelmParser& parser)
{
    if (!helm_version.empty() && helm_version != "V2.0") {
        parser.saveError(
            helm_version,
            "Only HELM and HELMV2.0 versions are currently supported");
    }

    if (helm_version.empty() && !extended_annotations.empty()) {
        parser.saveError(extended_annotations,
                         "HELM annotations are currently unsupported");
    }
}
} // namespace

void validate_parsed_info(const helm_info& parsed_info, HelmParser& parser)
{
    validate_helm_version(parsed_info.helm_version,
                          parsed_info.extended_annotations, parser);

    std::unordered_set<std::string_view> polymer_ids;
    std::for_each(parsed_info.polymers.begin(), parsed_info.polymers.end(),
                  [&](const auto& polymer) {
                      validate_polymer(polymer, polymer_ids, parser);
                  });
    if (parser.hasErrors()) {
        return;
    }

    std::for_each(parsed_info.connections.begin(),
                  parsed_info.connections.end(), [&](const auto& connection) {
                      validate_connection(connection, polymer_ids, parser);
                  });
    if (parser.hasErrors()) {
        return;
    }
    // if connection has no unsupported features check for missing residues
    // or incompatible rgroups
    std::for_each(parsed_info.connections.begin(),
                  parsed_info.connections.end(), [&](const auto& connection) {
                      validate_connection_residues(
                          connection, parsed_info.polymers, parser);
                  });
    if (parser.hasErrors()) {
        return;
    }

    std::unordered_set<std::string_view> group_ids;
    std::for_each(
        parsed_info.polymer_groups.begin(), parsed_info.polymer_groups.end(),
        [&](const auto& polymer_group) {
            validate_polymer_group_ids(polymer_group, group_ids, parser);
        });
    if (parser.hasErrors()) {
        return;
    }

    std::for_each(parsed_info.polymer_groups.begin(),
                  parsed_info.polymer_groups.end(),
                  [&](const auto& polymer_group) {
                      validate_polymer_group_items(polymer_group, polymer_ids,
                                                   group_ids, parser);
                  });

    validate_extended_annotations(parsed_info.extended_annotations, parser);
}
} // namespace helm
