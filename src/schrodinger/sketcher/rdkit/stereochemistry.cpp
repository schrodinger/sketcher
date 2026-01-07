#include "schrodinger/sketcher/rdkit/stereochemistry.h"

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/StereoGroup.h>

#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{
template <typename T> bool has_non_CIP_neighbor(const T& obj)
{
    auto stereo_info = RDKit::Chirality::detail::getStereoInfo(&obj);
    auto& mol = obj.getOwningMol();
    for (auto atom_idx : stereo_info.controllingAtoms) {
        // stereo bonds may have only 2 controlling atoms. The other 2
        // are NOATOM placeholders.
        if (atom_idx == RDKit::Chirality::StereoInfo::NOATOM) {
            continue;
        }
        auto controlling_atom = mol.getAtomWithIdx(atom_idx);
        if (RDKit::getAtomRLabel(controlling_atom) != 0 ||
            rdkit_extensions::is_attachment_point_dummy(*controlling_atom)) {
            return true;
        }
    }
    return false;
};
} // namespace

std::string get_atom_chirality_label(const RDKit::Atom& atom, bool strip_abs,
                                     bool show_unspecified)
{
    // SKETCH-1729: non-CIP ranked atoms have no priority; do not show any label
    if (has_non_CIP_neighbor(atom)) {
        return "";
    }

    auto possible_but_not_specified = [](const auto& atom) {
        int possible = 0;
        atom.getPropIfPresent(RDKit::common_properties::_ChiralityPossible,
                              possible);
        return !atom.hasProp(RDKit::common_properties::_CIPCode) && possible;
    };
    if (possible_but_not_specified(atom) && show_unspecified) {
        return "(?)";
    }

    std::string chirality;
    atom.getPropIfPresent<std::string>(RDKit::common_properties::atomNote,
                                       chirality);

    if (strip_abs &&
        chirality.find(rdkit_extensions::ABSOLUTE_STEREO_PREFIX) == 0) {
        // remove the rdkit_extensions::ABSOLUTE_STEREO_PREFIX
        chirality =
            chirality.substr(rdkit_extensions::ABSOLUTE_STEREO_PREFIX.size());
    }
    return chirality;
}

std::string get_bond_stereo_label(const RDKit::Bond& bond)
{
    std::string label;

    // SKETCH-2293: non-CIP ranked atoms have no priority; do not show any
    // label. has_non_CIP_neighbor () can't be called on non-stereo bonds
    // because getStereoInfo(bond) will throw if bond is not a stereo bond.
    if (bond.getPropIfPresent(RDKit::common_properties::bondNote, label) &&
        !label.empty() && has_non_CIP_neighbor(bond)) {
        return "";
    }

    return label;
}

/**
 * Read the enhanced stereo data, if any, for an atom.
 */
std::optional<EnhancedStereo>
get_enhanced_stereo_for_atom(const RDKit::Atom* const atom)
{
    auto& mol = atom->getOwningMol();
    for (const auto& stereo_group : mol.getStereoGroups()) {
        const auto& group_atoms = stereo_group.getAtoms();
        if (std::find(group_atoms.begin(), group_atoms.end(), atom) !=
            group_atoms.end()) {
            return EnhancedStereo(stereo_group.getGroupType(),
                                  stereo_group.getReadId());
        }
    }
    // this atom doesn't appear in any stereo groups
    return std::nullopt;
}

void set_enhanced_stereo_for_atom(RDKit::Atom* atom,
                                  const EnhancedStereo& enh_stereo)
{
    auto cur_stereo = get_enhanced_stereo_for_atom(atom);
    if (cur_stereo.has_value() && *cur_stereo == enh_stereo) {
        // the atom already has the desired stereo, so there's nothing to do
        return;
    }

    auto& mol = atom->getOwningMol();
    auto stereo_groups = mol.getStereoGroups();
    // remove the atom from whatever group it's currently in
    RDKit::removeAtomFromGroups(atom, stereo_groups);

    // check to see if there's already a stereo group with the desired
    // parameters
    std::vector<RDKit::Atom*> atoms;
    std::vector<RDKit::Bond*> bonds;
    bool is_abs = enh_stereo.first == RDKit::StereoGroupType::STEREO_ABSOLUTE;
    auto cur_group_it = stereo_groups.begin();
    for (; cur_group_it != stereo_groups.end(); ++cur_group_it) {
        auto& cur_group = *cur_group_it;
        if (cur_group.getGroupType() == enh_stereo.first &&
            (is_abs || enh_stereo.second == cur_group.getReadId())) {
            // this group matches the specified parameters
            atoms = cur_group.getAtoms();
            bonds = cur_group.getBonds();
            break;
        }
    }
    atoms.push_back(atom);
    auto group_id = is_abs ? 0 : enh_stereo.second;
    auto new_group =
        RDKit::StereoGroup(enh_stereo.first, atoms, bonds, group_id);
    new_group.setWriteId(group_id);
    if (cur_group_it == stereo_groups.end()) {
        // we didn't find any existing group with the required parameters, so we
        // add the new group at the end
        stereo_groups.push_back(new_group);
    } else {
        // We found a group in the for loop above. Because there's no public way
        // to add an atom to an existing stereo group instance, we instead
        // replace the group with a new one containing all of the same atoms and
        // bonds, plus the atom we want to add.
        *cur_group_it = new_group;
    }
    mol.setStereoGroups(stereo_groups);
}

} // namespace sketcher
} // namespace schrodinger
