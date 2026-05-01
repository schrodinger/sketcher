// -------------------------------------------------------------------------
// Implements polymer group APIs for HELM models
//
// Copyright Schrodinger LLC, All Rights Reserved.
// -------------------------------------------------------------------------
#include "schrodinger/rdkit_extensions/polymer_group.h"
#include "schrodinger/rdkit_extensions/helm.h"

#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ranges>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
std::vector<const RDKit::SubstanceGroup*>
get_polymer_group_substance_groups(const ::RDKit::ROMol& mol)
{
    std::vector<const RDKit::SubstanceGroup*> polymer_groups; // the result

    std::string sgroup_fieldname;
    std::string sgroup_type;
    std::string polymer_group_id;
    for (const auto& sgroup : ::RDKit::getSubstanceGroups(mol)) {
        if (!(sgroup.getPropIfPresent("TYPE", sgroup_type) &&
              sgroup.getPropIfPresent("FIELDNAME", sgroup_fieldname) &&
              sgroup.getPropIfPresent("ID", polymer_group_id))) {
            continue;
        }

        if (!(sgroup_type == "DAT" &&
              sgroup_fieldname.starts_with("POLYMER_GROUP") &&
              polymer_group_id.size() >= 2 && polymer_group_id[0] == 'G' &&
              std::ranges::all_of(polymer_group_id.substr(1),
                                  [](auto c) { return std::isdigit(c); }))) {
            continue;
        }

        polymer_groups.push_back(&sgroup);
    }

    return polymer_groups;
}

std::string get_polymer_group_helm_string(const ::RDKit::SubstanceGroup* sgroup)
{
    std::string id;
    std::string fieldname;
    std::vector<std::string> entities;
    std::vector<std::string> ratios;
    if (!sgroup->getPropIfPresent("ID", id)) {
        throw std::runtime_error(
            "Polymer group substance group missing required 'ID' property");
    }

    if (!sgroup->getPropIfPresent("FIELDNAME", fieldname)) {
        throw std::runtime_error(
            "Polymer group substance group missing required 'FIELDNAME' "
            "property");
    }

    if (!sgroup->getPropIfPresent(POLYMER_GROUP_ENTITIES, entities)) {
        throw std::runtime_error(fmt::format(
            "Polymer group substance group missing required '{}' property",
            POLYMER_GROUP_ENTITIES));
    }

    if (!sgroup->getPropIfPresent(POLYMER_GROUP_RATIOS, ratios)) {
        throw std::runtime_error(fmt::format(
            "Polymer group substance group missing required '{}' property",
            POLYMER_GROUP_RATIOS));
    }

    // Build HELM string from entities and ratios
    fmt::memory_buffer helm_string;
    fmt::format_to(std::back_inserter(helm_string), "{}({}", id, entities[0]);
    if (!ratios[0].empty()) {
        fmt::format_to(std::back_inserter(helm_string), ":{}", ratios[0]);
    }

    auto separator = (fieldname == POLYMER_GROUP_UNION ? '+' : ',');
    for (size_t i = 1; i < entities.size(); ++i) {
        fmt::format_to(std::back_inserter(helm_string), "{}{}", separator,
                       entities[i]);
        if (!ratios[i].empty()) {
            fmt::format_to(std::back_inserter(helm_string), ":{}", ratios[i]);
        }
    }

    fmt::format_to(std::back_inserter(helm_string), ")");
    return {helm_string.data(), helm_string.size()};
}

} // namespace

std::string get_polymer_groups_helm_string(const ::RDKit::ROMol& mol)
{
    // add this for backward compat for < 26-3
    if (auto sgroup = get_supplementary_info(mol); sgroup) {
        std::vector<std::string> supplementary_info;
        if (sgroup->getPropIfPresent("DATAFIELDS", supplementary_info) &&
            !supplementary_info.empty()) {
            return supplementary_info[0];
        }
        return {};
    }

    if (auto sgroups = get_polymer_group_substance_groups(mol);
        !sgroups.empty()) {
        return fmt::format(
            "{}", fmt::join(sgroups | std::views::transform([&](auto& sgroup) {
                                return get_polymer_group_helm_string(sgroup);
                            }),
                            "|"));
    }

    return {};
}

void add_polymer_group(RDKit::ROMol& mol, const std::string& name,
                       PolymerGroupType group_type,
                       const std::vector<std::string>& entities)
{
    // Delegate to overload with empty ratios
    std::vector<std::string> empty_ratios(entities.size(), "");
    add_polymer_group(mol, name, group_type, entities, empty_ratios);
}

void add_polymer_group(RDKit::ROMol& mol, const std::string& name,
                       PolymerGroupType group_type,
                       const std::vector<std::string>& entities,
                       const std::vector<std::string>& ratios)
{
    if (name.empty() || name[0] != 'G' || name.size() < 2 ||
        !std::ranges::all_of(name.substr(1),
                             [](char c) { return std::isdigit(c); })) {
        throw std::invalid_argument(
            "Polymer group name must start with 'G' followed by digits (e.g., "
            "'G1', 'G2')");
    }

    if (entities.size() < 2) {
        throw std::invalid_argument(
            "Polymer group must contain at least 2 entities");
    }

    if (ratios.size() != entities.size()) {
        throw std::invalid_argument(
            fmt::format("Ratios size ({}) must match entities size ({})",
                        ratios.size(), entities.size()));
    }

    std::unordered_map<std::string, std::vector<unsigned int>>
        available_entities;

    // add polymer groups
    for (auto sgroup : get_polymer_group_substance_groups(mol)) {
        auto polymer_group_id = sgroup->template getProp<std::string>("ID");
        available_entities[polymer_group_id] = sgroup->getAtoms();
    }

    for (auto atom : mol.atoms()) {
        if (atom->hasProp(REPETITION_DUMMY_ID)) { // query/repeated monomer
            continue;
        }

        auto polymer_id = get_polymer_id(atom);
        available_entities[polymer_id].push_back(atom->getIdx());
    }

    if (auto entry = available_entities.find(name);
        entry != available_entities.end()) {
        throw std::invalid_argument(
            fmt::format("Polymer group '{}' already exists", name));
    }

    if (auto entry = std::ranges::find_if(entities,
                                          [&](const auto& entity_id) {
                                              return available_entities.find(
                                                         entity_id) ==
                                                     available_entities.end();
                                          });
        entry != entities.end()) {
        throw std::invalid_argument(
            fmt::format("Entity '{}' not found in molecule", *entry));
    }

    // Collect unique atoms from all entities
    std::unordered_set<unsigned int> unique_atoms;
    for (const auto& entity_id : entities) {
        const auto& atoms = available_entities[entity_id];
        unique_atoms.insert(atoms.begin(), atoms.end());
    }

    // Convert to sorted vector for substance group
    std::vector<unsigned int> polymer_group_atoms(unique_atoms.begin(),
                                                  unique_atoms.end());
    std::ranges::sort(polymer_group_atoms);

    RDKit::SubstanceGroup sgroup{&mol, "DAT"};
    sgroup.setProp("FIELDNAME", group_type == PolymerGroupType::UNION
                                    ? POLYMER_GROUP_UNION
                                    : POLYMER_GROUP_EXCLUSIVE_LIST);
    sgroup.setProp("ID", name);

    // Store structured properties for efficient programmatic access
    sgroup.setProp(POLYMER_GROUP_ENTITIES, entities);
    sgroup.setProp(POLYMER_GROUP_RATIOS, ratios);

    sgroup.setAtoms(std::move(polymer_group_atoms));

    RDKit::addSubstanceGroup(mol, std::move(sgroup));
}

} // namespace rdkit_extensions
} // namespace schrodinger
