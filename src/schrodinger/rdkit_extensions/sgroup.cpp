#include "schrodinger/rdkit_extensions/sgroup.h"
#include <Geometry/point.h>
#include <GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string RDKIT_KEY_CONNECT = "CONNECT";
const std::string RDKIT_KEY_LABEL = "LABEL";
const std::string RDKIT_KEY_PARENT = "PARENT";
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

} // namespace rdkit_extensions
} // namespace schrodinger
