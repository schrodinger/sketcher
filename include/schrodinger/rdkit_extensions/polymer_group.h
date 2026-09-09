#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <optional>
#include <string>
#include <vector>

namespace RDKit
{
class ROMol;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

// Polymer group property names
const std::string POLYMER_GROUP{"POLYMER_GROUP"};
const std::string POLYMER_GROUP_HELM{"POLYMER_GROUP_HELM"};
const std::string POLYMER_GROUP_UNION{"POLYMER_GROUP_UNION"};
const std::string POLYMER_GROUP_EXCLUSIVE_LIST{"POLYMER_GROUP_EXCLUSIVE_LIST"};
const std::string POLYMER_GROUP_ENTITIES{"POLYMER_GROUP_ENTITIES"};
const std::string POLYMER_GROUP_RATIOS{"POLYMER_GROUP_RATIOS"};

/**
 * Type of polymer group grouping
 */
enum class PolymerGroupType {
    UNION,         // Polymers mixed together, separated by "+" in HELM
    EXCLUSIVE_LIST // Mutually exclusive alternatives, separated by "," in HELM
};

/**
 * Get the HELM v2.0 representation of all polymer groups in a molecule
 *
 * @param mol coarse grain molecule
 * @return the polymer groups as a HELM v2.0 string, or std::nullopt if none
 * exist
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::string
get_polymer_groups_helm_string(const ::RDKit::ROMol& mol);

/**
 * Add a polymer group to a monomeric molecule.
 *
 * Polymer groups represent hierarchical groupings of polymers or other polymer
 * groups, as defined in HELM v2.0 specification. They are stored as nested
 * substance groups to preserve hierarchy information during HELM round-tripping
 * and custom entity registration.
 *
 * Two types of polymer groups are supported:
 * - Union ("+") - polymers mixed together (e.g., antibody arms)
 * - Exclusive list (",") - mutually exclusive alternatives
 *
 * Groups can be nested by referencing other group IDs in the entities list.
 *
 * HELM v2.0 format examples:
 *   G1(PEPTIDE1+PEPTIDE2)              - Union of two polymers
 *   G1(PEPTIDE1,PEPTIDE2)              - Exclusive list
 *   G1(PEPTIDE1:2+PEPTIDE2:5)          - Union with ratios
 *   G1(PEPTIDE1+PEPTIDE2)|G2(G1+RNA1)  - Nested groups
 *
 * @param mol Monomeric molecule to add the polymer group to
 * @param name Polymer group identifier, must start with 'G' followed by digits
 *             (e.g., "G1", "G2"). Must be unique within the molecule.
 * @param group_type Type of polymer group (UNION for "+", EXCLUSIVE_LIST for
 * ",")
 * @param entities List of polymer IDs or group IDs to include in this group.
 *                 Must contain at least 2 entities.
 *
 * @throws std::invalid_argument if name doesn't start with 'G', if fewer than
 *         2 entities are provided, if the group name already exists, or if any
 *         entity is not found in the molecule
 */
RDKIT_EXTENSIONS_API void
add_polymer_group(RDKit::ROMol& mol, const std::string& name,
                  PolymerGroupType group_type,
                  const std::vector<std::string>& entities);

/**
 * Add a polymer group with ratios to a monomeric molecule.
 *
 * This overload allows specifying ratios for each entity in the polymer group.
 *
 * @param mol Monomeric molecule to add the polymer group to
 * @param name Polymer group identifier, must start with 'G' followed by digits
 *             (e.g., "G1", "G2"). Must be unique within the molecule.
 * @param group_type Type of polymer group (UNION for "+", EXCLUSIVE_LIST for
 * ",")
 * @param entities List of polymer IDs or group IDs to include in this group.
 *                 Must contain at least 2 entities.
 * @param ratios List of ratio strings for each entity. Use empty string for
 *               entities without ratios. Must be same length as entities.
 *               Examples: "2", "5", "0.1-1.3", ""
 *
 * @throws std::invalid_argument if name doesn't start with 'G', if fewer than
 *         2 entities are provided, if ratios size doesn't match entities size,
 *         if the group name already exists, or if any entity is not found
 */
RDKIT_EXTENSIONS_API void
add_polymer_group(RDKit::ROMol& mol, const std::string& name,
                  PolymerGroupType group_type,
                  const std::vector<std::string>& entities,
                  const std::vector<std::string>& ratios);

} // namespace rdkit_extensions
} // namespace schrodinger
