#include "schrodinger/rdkit_extensions/stereochemistry.h"

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/StereoGroup.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>

#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"

namespace schrodinger
{
namespace rdkit_extensions
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
            is_attachment_point_dummy(*controlling_atom)) {
            return true;
        }
    }
    return false;
};
} // namespace

UseModernStereoPerception::UseModernStereoPerception() :
    m_stereo_algo_state{RDKit::Chirality::getUseLegacyStereoPerception()}
{
    RDKit::Chirality::setUseLegacyStereoPerception(false);
}

UseModernStereoPerception::~UseModernStereoPerception()
{
    RDKit::Chirality::setUseLegacyStereoPerception(m_stereo_algo_state);
}

void assign_stereochemistry(RDKit::ROMol& mol)
{
    if (mol.getNumConformers() != 0) {
        auto conf = mol.getConformer();
        const auto& ps = conf.getPositions();

        if (std::any_of(ps.begin(), ps.end(),
                        [](auto p) { return p.length() > 1e-4; })) {
            // If we have all-zero coordinates, we don't care about
            // the 2D/3D flag (which is read, rather than checked)

            if (conf.is3D()) {
                // If we have 3D coordinates, calculate both chirality
                // and stereo bonds from them.
                RDKit::MolOps::assignStereochemistryFrom3D(mol);
                return;

            } else {
                // If we have 2D coordinates available, use them
                // to refresh stereo bond directions
                RDKit::MolOps::setDoubleBondNeighborDirections(mol, &conf);
            }
        }
    }

    // The general case, valid both no conformation and calculate
    // stereochemistry from parity, and stereo bonds from bond directions.
    bool cleanIt = false;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    RDKit::MolOps::assignStereochemistry(mol, cleanIt, force,
                                         flagPossibleStereoCenters);
}

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

    if (strip_abs && chirality.find(ABSOLUTE_STEREO_PREFIX) == 0) {
        // remove the ABSOLUTE_STEREO_PREFIX
        chirality = chirality.substr(ABSOLUTE_STEREO_PREFIX.size());
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

void wedgeMolBonds(RDKit::ROMol& mol, const RDKit::Conformer* conf)
{
    std::vector<RDKit::Bond*> attachment_dummy_bonds;
    for (auto atom : mol.atoms()) {
        if (is_attachment_point_dummy(*atom)) {
            auto bonds = mol.atomBonds(atom);
            auto bond = *bonds.begin();
            if (bond->getBondType() == RDKit::Bond::SINGLE) {
                attachment_dummy_bonds.push_back(bond);
                bond->setBondType(RDKit::Bond::OTHER);
            }
        }
    }

    RDKit::ClearSingleBondDirFlags(mol);
    RDKit::Chirality::clearMolBlockWedgingInfo(mol);

    try {
        // Temporarily silence RDKit's loggers
        RDLog::LogStateSetter silence_rdkit_logging;

        RDKit::WedgeMolBonds(mol, conf);
    } catch (const Invar::Invariant&) {
        // RDkit wasn't able to find a 'wedgeable' bond for some chiral atom
    }

    // Restore the dummies
    for (auto bond : attachment_dummy_bonds) {
        bond->setBondType(RDKit::Bond::SINGLE);
    }
}

/**
 * Read the enhanced stereo data, if any, for an atom.
 */
std::optional<rdkit_extensions::EnhancedStereo>
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

} // namespace rdkit_extensions
} // namespace schrodinger
