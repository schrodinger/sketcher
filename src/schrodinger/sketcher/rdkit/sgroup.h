#pragma once

#include "schrodinger/sketcher/definitions.h"
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <unordered_set>

namespace schrodinger
{
namespace sketcher
{

const float BRACKETS_LONG_SIDE = RDDepict::BOND_LEN * 0.9;

/**
 * Sets the coordinates of the corners of brackets to display for all S-groups
 * @param  mol the molecule for which all substance groups should be updated
 */
SKETCHER_API void update_s_group_brackets(RDKit::ROMol& mol);

/**
 * Removes a set of sgroups from a molecule
 * @param mol the molecule to remove the sgroups from
 * @param sgroups the set of sgroups to remove
 */
SKETCHER_API void remove_sgroups_from_molecule(
    RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups);

SKETCHER_API std::string get_sgroup_type(const RDKit::SubstanceGroup& sgroup);

SKETCHER_API std::string
get_sgroup_subtype(const RDKit::SubstanceGroup& sgroup);

SKETCHER_API std::string
get_repeat_pattern_label(const RDKit::SubstanceGroup& sgroup);

SKETCHER_API std::string get_polymer_label(const RDKit::SubstanceGroup& sgroup);

SKETCHER_API void set_sgroup_type(const RDKit::SubstanceGroup& sgroup,
                                  const std::string& value);

SKETCHER_API void set_sgroup_subtype(const RDKit::SubstanceGroup& sgroup,
                                     const std::string& value);

SKETCHER_API void set_repeat_pattern_label(const RDKit::SubstanceGroup& sgroup,
                                           const std::string& value);

SKETCHER_API void set_polymer_label(const RDKit::SubstanceGroup& sgroup,
                                    const std::string& value);

/**
 * If the molecule contains an S-group that covers exactly the specified atoms,
 * then return that S-group.  Otherwise, return nullptr.
 */
SKETCHER_API const RDKit::SubstanceGroup*
get_existing_sgroup_for_atoms(std::unordered_set<const RDKit::Atom*> atoms,
                              const RDKit::ROMol& mol);

/**
 * Determine whether the specified atoms can form an S-group.  A group of atoms
 * can form an S-group if the following two conditions are met:
 * - the atoms are contiguous
 * - there are exactly two bonds between the specified atoms and all other atoms
 *   in the molecule
 */
SKETCHER_API bool
can_atoms_form_sgroup(std::unordered_set<const RDKit::Atom*> specified_atoms,
                      const RDKit::ROMol& mol);

/**
 * @return the atoms contained in the given S-group
 */
SKETCHER_API std::unordered_set<const RDKit::Atom*>
get_sgroup_atoms(const RDKit::SubstanceGroup* const s_group,
                 const RDKit::ROMol& mol);

/**
 * Find the bond indices of the two bonds that cross between the specified
 * atoms and all other atoms in the molecule.  This methods assumes that the
 * specified atoms are contiguous and that there are exactly two such bonds.  In
 * other words, can_atoms_form_s_group() should return true for these atoms.
 * This method *may* throw otherwise, but is not guaranteed to.
 */
SKETCHER_API std::vector<unsigned int>
get_bonds_for_sgroup_atoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                           const RDKit::ROMol& mol);

/**
 * @return all bonds between the atoms contained in the S-group
 */
SKETCHER_API std::unordered_set<const RDKit::Bond*>
get_bonds_within_sgroup(const RDKit::SubstanceGroup& s_group);

} // namespace sketcher
} // namespace schrodinger
