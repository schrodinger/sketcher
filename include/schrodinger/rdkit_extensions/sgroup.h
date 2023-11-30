#pragma once

#include "schrodinger/rdkit_extensions/definitions.h" // RDKIT_EXTENSIONS_API
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <unordered_set>

namespace schrodinger
{
namespace rdkit_extensions
{

const float BRACKETS_LONG_SIDE = RDDepict::BOND_LEN * 0.9;

/**
 * add the coordinates of the corners of brackets to display for an
 * S-group
 * @param  sgroup the substance group to add brackets to
 */
RDKIT_EXTENSIONS_API void
update_s_group_brackets(RDKit::SubstanceGroup& sgroup);

/**
 * remove a set of sgroups from a molecule
 * @param mol the molecule to remove the sgroups from
 * @param sgroups the set of sgroups to remove
 */
RDKIT_EXTENSIONS_API void remove_sgroups_from_molecule(
    RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups);

RDKIT_EXTENSIONS_API std::string
get_sgroup_type(const RDKit::SubstanceGroup& sgroup);

RDKIT_EXTENSIONS_API std::string
get_sgroup_subtype(const RDKit::SubstanceGroup& sgroup);

RDKIT_EXTENSIONS_API std::string
get_repeat_pattern_label(const RDKit::SubstanceGroup& sgroup);

RDKIT_EXTENSIONS_API std::string
get_polymer_label(const RDKit::SubstanceGroup& sgroup);

} // namespace rdkit_extensions
} // namespace schrodinger
