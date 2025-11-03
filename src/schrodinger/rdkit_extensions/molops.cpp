/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: miscellaneous mol operations
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/molops.h"

#include <fmt/format.h>
#include <rdkit/GraphMol/ChemTransforms/ChemTransforms.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
struct [[nodiscard]] selected_atom_info {
    std::vector<bool> selected_atoms;
    std::vector<bool> selected_bonds;
    std::unordered_map<unsigned int, unsigned int> atom_mapping;
    std::unordered_map<unsigned int, unsigned int> bond_mapping;
};

} // namespace

static void validate_monomeric_mol_input(const RDKit::ROMol&,
                                         const std::string_view&);

[[nodiscard]] static std::unordered_map<char, int>
get_max_polymer_suffixes(const std::vector<std::string>&);

[[nodiscard]] static std::unordered_map<std::string, std::string>
get_renamed_polymer_ids(const std::vector<std::string>&,
                        std::unordered_map<char, int>&);

void apply_sanitization(RDKit::RWMol& mol, Sanitization sanitization)
{
    RDLog::LogStateSetter silence_rdkit_logging;

    unsigned failed_op{0};
    int ops = RDKit::MolOps::SANITIZE_ALL;
    // Partial sanitization runs all operations except:
    // - cleanup: to avoid altering the chemistry of the input,
    // - properties: so we don't freak out about invalid valences,
    // - aromaticity and kekulization: to preserve the state the input comes in,
    // - chirality cleanup: to avoid changing stereo of the input.
    // We do it this way so we get any new ops that might be added to
    // SANITIZE_ALL.
    if (sanitization == Sanitization::PARTIAL) {
        ops &= ~(RDKit::MolOps::SANITIZE_CLEANUP |
                 RDKit::MolOps::SANITIZE_PROPERTIES |
                 RDKit::MolOps::SANITIZE_KEKULIZE |
                 RDKit::MolOps::SANITIZE_SETAROMATICITY |
                 RDKit::MolOps::SANITIZE_CLEANUPCHIRALITY);
    }
    RDKit::MolOps::sanitizeMol(mol, failed_op, ops);
    // Regardless of sanitization level, ensure property cache is updated
    mol.updatePropertyCache(false);
}

void addHs(RDKit::RWMol& mol, std::vector<unsigned> atom_ids)
{
    // If atom_ids is empty, add Hs to all atoms
    auto only_on_atoms = atom_ids.empty() ? nullptr : &atom_ids;

    auto initial_num_atoms = mol.getNumAtoms();
    bool explicit_only = false;
    bool add_coords = false;
    RDKit::MolOps::addHs(mol, explicit_only, add_coords, only_on_atoms);

    // Explicitly add 2D coordinates to the new hydrogens; ids are guaranteed to
    // be from the previous number of atoms to the current number of atoms
    std::vector<unsigned int> frozen_ids;
    for (unsigned int idx = 0; idx < initial_num_atoms; ++idx) {
        frozen_ids.push_back(idx);
    }
    rdkit_extensions::compute2DCoords(mol, frozen_ids);
}

void removeHs(RDKit::RWMol& rdk_mol)
{
    RDKit::MolOps::RemoveHsParameters ps;

    // We always remove H on queries; for sketcher import, all atoms are created
    // as QueryAtoms as that they might be changed into queries later on; for
    // conversion from 3D Structure, there is no way to create queries.
    ps.removeWithQuery = true;

    // preserve hydrides
    ps.removeHydrides = false;

    // Disable displaying warnings
    ps.showWarnings = false;

    bool sanitize = false;
    RDKit::MolOps::removeHs(rdk_mol, ps, sanitize);

    rdk_mol.updatePropertyCache(false);
}

void removeHs(RDKit::RWMol& rdk_mol, std::vector<unsigned> atom_ids)
{
    if (atom_ids.empty()) {
        return;
    }

    // Augment atom ids with the ids of the hydrogens attached to them.
    // We don't care about duplicates, because we will sort ant uniquify
    // the list later on.
    const auto num_atom_ids = atom_ids.size();
    for (unsigned i = 0; i < num_atom_ids; ++i) {
        const auto atom = rdk_mol.getAtomWithIdx(atom_ids[i]);
        for (auto nbr : rdk_mol.atomNeighbors(atom)) {
            if (nbr->getAtomicNum() == 1 && nbr->getIsotope() == 0) {
                atom_ids.push_back(nbr->getIdx());
            }
        }
    }

    constexpr int h_protection_mark = 1000;

    // Sort and uniquify the list of atom ids.
    std::sort(atom_ids.begin(), atom_ids.end());
    atom_ids.erase(std::unique(atom_ids.begin(), atom_ids.end()),
                   atom_ids.end());

    auto remove_id_itr = atom_ids.begin();
    std::vector<RDKit::Atom*> protected_atoms;
    const auto input_num_atoms = rdk_mol.getNumAtoms();
    for (unsigned i = 0; i < input_num_atoms; ++i) {
        if (remove_id_itr != atom_ids.end() && i == *remove_id_itr) {
            ++remove_id_itr;
            continue;
        }
        auto atom = rdk_mol.getAtomWithIdx(i);
        if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
            atom->setIsotope(atom->getIsotope() + h_protection_mark);
            protected_atoms.push_back(atom);
        }
    }

    rdkit_extensions::removeHs(rdk_mol);

    for (auto atom : protected_atoms) {
        atom->setIsotope(atom->getIsotope() - h_protection_mark);
    }
}

static void copy_selected_atoms_and_bonds(::RDKit::RWMol& extracted_mol,
                                          const RDKit::ROMol& reference_mol,
                                          selected_atom_info& selection_info)
{
    auto& [selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
        selection_info;
    for (const auto& ref_atom : reference_mol.atoms()) {
        if (!selected_atoms[ref_atom->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Atom> extracted_atom{ref_atom->copy()};

        static constexpr bool updateLabel = true;
        static constexpr bool takeOwnership = true;
        atom_mapping[ref_atom->getIdx()] = extracted_mol.addAtom(
            extracted_atom.release(), updateLabel, takeOwnership);
    }

    for (const auto& ref_bond : reference_mol.bonds()) {
        if (!selected_bonds[ref_bond->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Bond> extracted_bond{ref_bond->copy()};
        extracted_bond->setBeginAtomIdx(
            atom_mapping[ref_bond->getBeginAtomIdx()]);
        extracted_bond->setEndAtomIdx(atom_mapping[ref_bond->getEndAtomIdx()]);

        static constexpr bool takeOwnership = true;
        auto num_bonds =
            extracted_mol.addBond(extracted_bond.release(), takeOwnership);
        bond_mapping[ref_bond->getIdx()] = num_bonds - 1;
    }
}

[[nodiscard]] static bool
is_selected_sgroup(const ::RDKit::SubstanceGroup& sgroup,
                   const selected_atom_info& selection_info)
{
    auto is_selected_component = [](auto& indices, auto& selection_test) {
        return indices.empty() ||
               std::all_of(indices.begin(), indices.end(), selection_test);
    };

    // clang-format off
    auto atom_test = [&](int idx) { return selection_info.selected_atoms[idx]; };
    auto bond_test = [&](int idx) { return selection_info.selected_bonds[idx]; };
    return is_selected_component(sgroup.getAtoms(), atom_test) &&
           is_selected_component(sgroup.getBonds(), bond_test) &&
           is_selected_component(sgroup.getParentAtoms(), atom_test);
    // clang-format on
}

static void
copy_selected_substance_groups(::RDKit::RWMol& extracted_mol,
                               const RDKit::ROMol& reference_mol,
                               const selected_atom_info& selection_info)
{
    auto update_indices = [](auto& sgroup, auto getter, auto setter,
                             auto& mapping) {
        auto indices = getter(sgroup);
        std::for_each(indices.begin(), indices.end(),
                      [&](auto& idx) { idx = mapping.at(idx); });
        setter(sgroup, std::move(indices));
    };

    const auto& [selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
        selection_info;
    for (const auto& sgroup : ::RDKit::getSubstanceGroups(reference_mol)) {
        if (!is_selected_sgroup(sgroup, selection_info)) {
            continue;
        }

        ::RDKit::SubstanceGroup extracted_sgroup(sgroup);
        extracted_sgroup.setOwningMol(&extracted_mol);

        update_indices(
            extracted_sgroup, std::mem_fn(&::RDKit::SubstanceGroup::getAtoms),
            std::mem_fn(&::RDKit::SubstanceGroup::setAtoms), atom_mapping);
        update_indices(extracted_sgroup,
                       std::mem_fn(&::RDKit::SubstanceGroup::getParentAtoms),
                       std::mem_fn(&::RDKit::SubstanceGroup::setParentAtoms),
                       atom_mapping);
        update_indices(
            extracted_sgroup, std::mem_fn(&::RDKit::SubstanceGroup::getBonds),
            std::mem_fn(&::RDKit::SubstanceGroup::setBonds), bond_mapping);

        ::RDKit::addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
    }
}

static void
copy_selected_stereo_groups(::RDKit::RWMol& extracted_mol,
                            const RDKit::ROMol& reference_mol,
                            const selected_atom_info& selection_info)
{
    auto is_selected_component = [](auto& objects, auto& selected_indices) {
        return objects.empty() ||
               std::all_of(objects.begin(), objects.end(), [&](auto& object) {
                   return selected_indices[object->getIdx()];
               });
    };

    auto is_selected_stereo_group = [&](const auto& stereo_group) {
        return is_selected_component(stereo_group.getAtoms(),
                                     selection_info.selected_atoms) &&
               is_selected_component(stereo_group.getBonds(),
                                     selection_info.selected_bonds);
    };

    std::vector<::RDKit::Atom*> extracted_atoms(extracted_mol.getNumAtoms());
    for (const auto& atom : extracted_mol.atoms()) {
        extracted_atoms[atom->getIdx()] = atom;
    }

    std::vector<::RDKit::Bond*> extracted_bonds(extracted_mol.getNumBonds());
    for (const auto& bond : extracted_mol.bonds()) {
        extracted_bonds[bond->getIdx()] = bond;
    }

    const auto& [selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
        selection_info;
    std::vector<::RDKit::StereoGroup> extracted_stereo_groups;
    for (const auto& stereo_group : reference_mol.getStereoGroups()) {
        if (!is_selected_stereo_group(stereo_group)) {
            continue;
        }

        std::vector<::RDKit::Atom*> atoms;
        for (const auto& atom : stereo_group.getAtoms()) {
            atoms.push_back(extracted_atoms[atom_mapping.at(atom->getIdx())]);
        }

        std::vector<::RDKit::Bond*> bonds;
        for (const auto& bond : stereo_group.getBonds()) {
            bonds.push_back(extracted_bonds[bond_mapping.at(bond->getIdx())]);
        }

        extracted_stereo_groups.push_back({stereo_group.getGroupType(),
                                           std::move(atoms), std::move(bonds),
                                           stereo_group.getReadId()});
        extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
    }

    extracted_mol.setStereoGroups(std::move(extracted_stereo_groups));
}

boost::shared_ptr<RDKit::RWMol>
ExtractMolFragment(const RDKit::ROMol& mol,
                   const std::vector<unsigned int>& atom_ids, bool sanitize)
{

    const auto num_atoms = mol.getNumAtoms();
    selected_atom_info selection_info{std::vector<bool>(num_atoms),
                                      std::vector<bool>(mol.getNumBonds()),
                                      {},
                                      {}};
    auto& [selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
        selection_info;
    for (const auto& atom_idx : atom_ids) {
        if (atom_idx < num_atoms) {
            selected_atoms[atom_idx] = true;
        }
    }
    for (const auto& bond : mol.bonds()) {
        if (selected_atoms[bond->getBeginAtomIdx()] &&
            selected_atoms[bond->getEndAtomIdx()]) {
            selected_bonds[bond->getIdx()] = true;
        }
    }

    auto extracted_mol = std::make_unique<::RDKit::RWMol>();
    copy_selected_atoms_and_bonds(*extracted_mol, mol, selection_info);
    copy_selected_substance_groups(*extracted_mol, mol, selection_info);
    copy_selected_stereo_groups(*extracted_mol, mol, selection_info);
    if (sanitize) {
        ::RDKit::MolOps::sanitizeMol(*extracted_mol);
    }

    // NOTE: Bookmarks are currently not copied
    return extracted_mol;
}

boost::shared_ptr<RDKit::ROMol> CombineMonomericMols(const RDKit::ROMol& mol1,
                                                     const RDKit::ROMol& mol2)
{
    // make sure there are no unsupported features
    validate_monomeric_mol_input(mol1, "mol1");
    validate_monomeric_mol_input(mol2, "mol2");

    std::unique_ptr<RDKit::ROMol> result;

    auto mol1_polymer_ids = get_polymer_ids(mol1);
    auto mol2_polymer_ids = get_polymer_ids(mol2);
    // We need to de-duplicate the polymer ids
    if (std::ranges::find_first_of(mol1_polymer_ids, mol2_polymer_ids) !=
        mol1_polymer_ids.end()) {
        // Sort second set of ids so we see POLYMER1 before POLYMER2
        std::ranges::sort(mol2_polymer_ids);

        // rename polymer ids to de-duplicate them
        auto max_polymer_suffixes = get_max_polymer_suffixes(mol1_polymer_ids);
        auto renamed_polymer_ids =
            get_renamed_polymer_ids(mol2_polymer_ids, max_polymer_suffixes);

        // create mutable copy
        RDKit::RWMol mol2_copy(mol2);

        // Atom-level updates:
        //      * REPETITION_DUMMY_ID
        //      * atom chain id
        std::string val;
        for (auto& atom : mol2_copy.atoms()) {
            if (atom->getPropIfPresent(REPETITION_DUMMY_ID, val)) {
                atom->setProp(REPETITION_DUMMY_ID, renamed_polymer_ids[val]);
            }

            if (auto res_info = static_cast<RDKit::AtomPDBResidueInfo*>(
                    atom->getMonomerInfo());
                res_info) {
                res_info->setChainId(
                    renamed_polymer_ids[res_info->getChainId()]);
            }
        }

        // Sgroup-level updates:
        //      * polymer annotation SubstanceGroup ID
        for (auto& sgroup : RDKit::getSubstanceGroups(mol2_copy)) {
            if (sgroup.getPropIfPresent("TYPE", val) && val == "COP" &&
                sgroup.getPropIfPresent("ID", val)) {
                sgroup.setProp("ID", renamed_polymer_ids[val]);
            }
        }

        result.reset(RDKit::combineMols(mol1, mol2_copy));
    } else {
        result.reset(RDKit::combineMols(mol1, mol2));
    }

    return result;
}

static void validate_monomeric_mol_input(const RDKit::ROMol& mol,
                                         const std::string_view& name)
{
    if (!isMonomeric(mol)) {
        auto msg = fmt::format("'{}' must be a monomeric mol", name);
        throw std::invalid_argument(msg);
    }

    if (mol.hasProp(EXTENDED_ANNOTATIONS)) {
        auto msg = fmt::format(
            "'{}' has unsupported feature: EXTENDED_ANNOTATIONS", name);
        throw std::invalid_argument(msg);
    }

    if (const auto sgroup = get_supplementary_info(mol); sgroup) {
        std::vector<std::string> datafields;
        if (!sgroup->getPropIfPresent("DATAFIELDS", datafields)) {
            return;
        }

        if (!datafields[0].empty()) {
            auto msg = fmt::format(
                "'{}' has unsupported feature: POLYMER_GROUPS", name);
            throw std::invalid_argument(msg);
        }

        if (!datafields[1].empty()) {
            auto msg = fmt::format(
                "'{}' has unsupported feature: EXTENDED_ANNOTATIONS", name);
            throw std::invalid_argument(msg);
        }
    }
}

[[nodiscard]] static std::unordered_map<char, int>
get_max_polymer_suffixes(const std::vector<std::string>& polymer_ids)
{
    std::unordered_map<char, int> max_polymer_suffixes;
    for (auto& polymer_id : polymer_ids) {
        auto pos = polymer_id.find_first_of("0123456789");
        auto suffix = std::stoi(polymer_id.substr(pos));
        if (suffix > max_polymer_suffixes[polymer_id[0]]) {
            max_polymer_suffixes[polymer_id[0]] = suffix;
        }
    }

    return max_polymer_suffixes;
}

[[nodiscard]] static std::unordered_map<std::string, std::string>
get_renamed_polymer_ids(const std::vector<std::string>& polymer_ids,
                        std::unordered_map<char, int>& max_polymer_suffixes)
{
    std::unordered_map<std::string, std::string> renamed_polymer_ids;
    for (const auto& polymer_id : polymer_ids) {
        // assume we don't need to de-deplicate
        if (max_polymer_suffixes[polymer_id[0]] == 0) {
            renamed_polymer_ids[polymer_id] = polymer_id;
        }
        // rename polymer id by updating suffix
        else {
            auto pos = polymer_id.find_first_of("0123456789");
            auto prefix = polymer_id.substr(0, pos);
            auto renamed_polymer_id = fmt::format(
                "{}{}", prefix, ++max_polymer_suffixes[polymer_id[0]]);
            renamed_polymer_ids[polymer_id] = renamed_polymer_id;
        }
    }
    return renamed_polymer_ids;
}

boost::shared_ptr<RDKit::ROMol> CombineMols(const RDKit::ROMol& mol1,
                                            const RDKit::ROMol& mol2)
{
    if (!isMonomeric(mol1) && !isMonomeric(mol2)) {
        return std::unique_ptr<RDKit::ROMol>{RDKit::combineMols(mol1, mol2)};
    } else if (!isMonomeric(mol1)) {
        auto monomermol1 = toMonomeric(mol1);
        return CombineMonomericMols(*monomermol1, mol2);
    } else if (!isMonomeric(mol2)) {
        auto monomermol2 = toMonomeric(mol2);
        return CombineMonomericMols(mol1, *monomermol2);
    } else {
        return CombineMonomericMols(mol1, mol2);
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
