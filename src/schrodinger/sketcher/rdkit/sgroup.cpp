#include "schrodinger/sketcher/rdkit/sgroup.h"
#include <Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

namespace schrodinger
{
namespace sketcher
{

const std::string RDKIT_KEY_CONNECT = "CONNECT";
const std::string RDKIT_KEY_LABEL = "LABEL";
const std::string RDKIT_KEY_SUBTYPE = "SUBTYPE";
const std::string RDKIT_KEY_TYPE = "TYPE";

/**
 * @param sgroup the SGroup for which to update bracket coordinates
 */
static void update_bracket_coordinates(RDKit::SubstanceGroup& sgroup)
{
    if (!sgroup.hasOwningMol()) {
        throw std::runtime_error("S-group does not belong to a molecule");
    }
    sgroup.clearBrackets();
    auto molecule = sgroup.getOwningMol();
    // get all connection bonds, i.e. bonds that connect one atom of the SGroup
    // to one that is outside of it
    auto connection_bonds = sgroup.getBonds();
    // only add brackets if there's exactly two connection bonds
    if (connection_bonds.size() != 2) {
        return;
    }
    for (auto idx : connection_bonds) {
        auto connection_bond = molecule.getBondWithIdx(idx);
        auto atom1 = connection_bond->getBeginAtomIdx();
        auto atom2 = connection_bond->getEndAtomIdx();
        auto pos1 = molecule.getConformer().getAtomPos(atom1);
        auto pos2 = molecule.getConformer().getAtomPos(atom2);

        auto mid_point = (pos1 + pos2) * 0.5;
        auto bond_dir = pos1 - pos2;
        auto normal = RDGeom::Point3D(-bond_dir.y, bond_dir.x, 0.0);
        normal.normalize();
        normal *= BRACKETS_LONG_SIDE * 0.5;
        sgroup.addBracket(
            {mid_point + normal, mid_point - normal, RDGeom::Point3D(0, 0, 0)});
    }
}

void update_s_group_brackets(RDKit::ROMol& mol)
{
    for (auto& sgroup : getSubstanceGroups(mol)) {
        update_bracket_coordinates(sgroup);
    }
}

void remove_sgroups_from_molecule(
    RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups)
{
    // RDKit doesn't provide a way to remove sgroups from a molecule, so we have
    // to use this workaround
    auto& mol_sgroups = getSubstanceGroups(mol);
    std::vector<RDKit::SubstanceGroup> new_sgroups;
    for (auto& sgroup : mol_sgroups) {
        if (std::find(sgroups.begin(), sgroups.end(), &sgroup) ==
            sgroups.end()) {
            new_sgroups.push_back(sgroup);
        }
    }
    mol_sgroups = std::move(new_sgroups);
}

static std::string get_string_property(const RDKit::SubstanceGroup& sgroup,
                                       const std::string& key)
{
    std::string value;
    sgroup.getPropIfPresent(key, value);
    return value;
}

std::string get_sgroup_type(const RDKit::SubstanceGroup& sgroup)
{
    return get_string_property(sgroup, RDKIT_KEY_TYPE);
}

std::string get_sgroup_subtype(const RDKit::SubstanceGroup& sgroup)
{
    return get_string_property(sgroup, RDKIT_KEY_SUBTYPE);
}

std::string get_repeat_pattern_label(const RDKit::SubstanceGroup& sgroup)
{
    return get_string_property(sgroup, RDKIT_KEY_CONNECT);
}

std::string get_polymer_label(const RDKit::SubstanceGroup& sgroup)
{
    return get_string_property(sgroup, RDKIT_KEY_LABEL);
}

static void set_string_property(const RDKit::SubstanceGroup& sgroup,
                                const std::string& key,
                                const std::string& value)
{
    if (value.empty()) {
        sgroup.clearProp(key);
    } else {
        sgroup.setProp(key, value);
    }
}

void set_sgroup_type(const RDKit::SubstanceGroup& sgroup,
                     const std::string& value)
{
    set_string_property(sgroup, RDKIT_KEY_TYPE, value);
}

void set_sgroup_subtype(const RDKit::SubstanceGroup& sgroup,
                        const std::string& value)
{
    set_string_property(sgroup, RDKIT_KEY_SUBTYPE, value);
}

void set_repeat_pattern_label(const RDKit::SubstanceGroup& sgroup,
                              const std::string& value)
{
    set_string_property(sgroup, RDKIT_KEY_CONNECT, value);
}

void set_polymer_label(const RDKit::SubstanceGroup& sgroup,
                       const std::string& value)
{
    set_string_property(sgroup, RDKIT_KEY_LABEL, value);
}

const RDKit::SubstanceGroup*
get_existing_sgroup_for_atoms(std::unordered_set<const RDKit::Atom*> atoms,
                              const RDKit::ROMol& mol)
{
    for (auto& s_group : getSubstanceGroups(mol)) {
        auto s_group_atoms = get_sgroup_atoms(&s_group, mol);
        if (s_group_atoms == atoms) {
            return &s_group;
        }
    }
    return nullptr;
}

bool can_atoms_form_sgroup(
    std::unordered_set<const RDKit::Atom*> specified_atoms,
    const RDKit::ROMol& mol)
{
    if (specified_atoms.empty()) {
        return false;
    }
    std::unordered_set<const RDKit::Atom*> atoms_to_visit,
        visited_specified_atoms, visited_unspecified_atoms;
    const RDKit::Atom* starting_atom = *specified_atoms.begin();
    atoms_to_visit.insert(starting_atom);
    while (!atoms_to_visit.empty()) {
        auto* cur_atom = *atoms_to_visit.begin();
        atoms_to_visit.erase(cur_atom);
        if (specified_atoms.count(cur_atom)) {
            // if the atom is specified, then we want to visit all of its
            // neighbors
            visited_specified_atoms.insert(cur_atom);
            for (auto* cur_neighbor : mol.atomNeighbors(cur_atom)) {
                if (!(atoms_to_visit.count(cur_neighbor) ||
                      visited_specified_atoms.count(cur_neighbor) ||
                      visited_unspecified_atoms.count(cur_neighbor))) {
                    atoms_to_visit.insert(cur_neighbor);
                }
            }
        } else {
            // if the atom isn't specified, then we've found a bond between
            // specified and unspecified atoms
            visited_unspecified_atoms.insert(cur_atom);
            if (visited_unspecified_atoms.size() > 2) {
                // there are more than two bonds that cross between specified
                // and unspecified atoms, so these atoms can't form an S-group
                return false;
            }
        }
    }
    // make sure that we've visited all specified atoms and exactly two
    // unspecified atoms
    return (visited_specified_atoms.size() == specified_atoms.size() &&
            visited_unspecified_atoms.size() == 2);
}

std::unordered_set<const RDKit::Atom*>
get_sgroup_atoms(const RDKit::SubstanceGroup* const s_group,
                 const RDKit::ROMol& mol)
{
    const auto& s_group_atom_idxs = s_group->getAtoms();
    std::unordered_set<const RDKit::Atom*> s_group_atoms;
    std::transform(s_group_atom_idxs.begin(), s_group_atom_idxs.end(),
                   std::inserter(s_group_atoms, s_group_atoms.begin()),
                   [&mol](auto idx) { return mol.getAtomWithIdx(idx); });
    return s_group_atoms;
}

std::vector<unsigned int>
get_bonds_for_sgroup_atoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                           const RDKit::ROMol& mol)
{
    std::vector<unsigned int> bond_idxs;
    for (auto* cur_atom : atoms) {
        for (auto* cur_neighbor : mol.atomNeighbors(cur_atom)) {
            if (!atoms.count(cur_neighbor)) {
                auto* bond = mol.getBondBetweenAtoms(cur_atom->getIdx(),
                                                     cur_neighbor->getIdx());
                bond_idxs.push_back(bond->getIdx());
                if (bond_idxs.size() == 2) {
                    return bond_idxs;
                }
            }
        }
    }
    throw std::runtime_error("Could not find two S-group bonds");
}

std::unordered_set<const RDKit::Bond*>
get_bonds_within_sgroup(const RDKit::SubstanceGroup& s_group)
{
    std::unordered_set<const RDKit::Bond*> bonds;
    std::unordered_set<int> s_group_atom_idxs(s_group.getAtoms().begin(),
                                              s_group.getAtoms().end());
    auto& mol = s_group.getOwningMol();
    for (auto atom_idx : s_group.getAtoms()) {
        auto* atom = mol.getAtomWithIdx(atom_idx);
        for (auto* cur_bond : mol.atomBonds(atom)) {
            unsigned int other_atom_idx = cur_bond->getOtherAtomIdx(atom_idx);
            if (s_group_atom_idxs.count(other_atom_idx)) {
                bonds.insert(cur_bond);
            }
        }
    }
    return bonds;
}

} // namespace sketcher
} // namespace schrodinger
