// -------------------------------------------------------------------------
// Implements schrodinger::rdkit_extensions:: helper apis for HELM models
//
//
// Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/helm.h"

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/container/dynamic_bitset_on_bits_wrapper.hpp"
#include "schrodinger/rdkit_extensions/molops.h"

using ::RDKit::RWMol;

namespace schrodinger
{
namespace rdkit_extensions
{

[[nodiscard]] bool is_coarse_grain_mol(const RDKit::ROMol& mol)
{
    return mol.hasProp(HELM_MODEL);
}

static void check_if_input_is_atomistic(const RDKit::ROMol& mol)
{
    if (!is_coarse_grain_mol(mol)) {
        throw std::invalid_argument("Atomistic mols are currently unsupported");
    }
}

[[nodiscard]] static const RDKit::SubstanceGroup*
find_polymer_sgroup(const RDKit::ROMol& mol, const std::string_view& target_id)
{
    const auto& sgroups = ::RDKit::getSubstanceGroups(mol);

    std::string polymer_type;
    std::string polymer_id;
    auto sgroup =
        std::find_if(sgroups.begin(), sgroups.end(), [&](auto& sgroup) {
            return sgroup.getPropIfPresent("TYPE", polymer_type) &&
                   sgroup.getPropIfPresent("ID", polymer_id) &&
                   polymer_type == "COP" && polymer_id == target_id;
        });

    return sgroup == sgroups.end() ? nullptr : &(*sgroup);
}

[[nodiscard]] std::vector<unsigned int>
get_atoms_in_polymer_chain(const RDKit::ROMol& mol, std::string_view polymer_id)
{
    check_if_input_is_atomistic(mol);

    auto polymer = find_polymer_sgroup(mol, polymer_id);
    return polymer ? polymer->getAtoms() : std::vector<unsigned int>();
}

[[nodiscard]] std::vector<unsigned int>
get_atoms_in_polymer_chains(const RDKit::ROMol& mol,
                            const std::vector<std::string_view>& polymer_ids)
{
    check_if_input_is_atomistic(mol);

    boost::dynamic_bitset<> selected_atoms(mol.getNumAtoms());
    for (auto& polymer_id : polymer_ids) {
        auto polymer = find_polymer_sgroup(mol, polymer_id);
        if (!polymer) {
            continue;
        }

        for (auto& atom_idx : polymer->getAtoms()) {
            selected_atoms.set(atom_idx);
        }
    }

    schrodinger::dynamic_bitset_on_bits_wrapper on_bits(selected_atoms);
    return {on_bits.begin(), on_bits.end()};
}

static void
check_for_polymer_groups_and_extended_annotations(const ::RDKit::ROMol& mol)
{
    std::string type;
    std::string fieldname;
    for (auto& sgroup : ::RDKit::getSubstanceGroups(mol)) {
        static const std::string TYPE{"TYPE"};
        static const std::string FIELDNAME{"FIELDNAME"};
        if (!(sgroup.getPropIfPresent(TYPE, type) &&
              sgroup.getPropIfPresent(FIELDNAME, fieldname) &&
              fieldname == SUPPLEMENTARY_INFORMATION)) {
            continue;
        }

        std::vector<std::string> datafields;
        static const std::string DATAFIELDS{"DATAFIELDS"};
        if (!sgroup.getPropIfPresent(DATAFIELDS, datafields) ||
            datafields.size() < 2u) {
            continue;
        }

        if (!datafields[0].empty()) {
            throw std::invalid_argument(
                "Polymer groups are currently unsupported");
        }

        if (!datafields[1].empty()) {
            throw std::invalid_argument(
                "Extended annotations are currently unsupported");
        }
        return;
    }
}

[[nodiscard]] boost::shared_ptr<RDKit::RWMol>
extract_helm_polymers(const RDKit::ROMol& mol,
                      const std::vector<std::string_view>& polymer_ids)
{
    check_if_input_is_atomistic(mol);
    check_for_polymer_groups_and_extended_annotations(mol);

    auto selected_atom_indices = get_atoms_in_polymer_chains(mol, polymer_ids);
    constexpr bool sanitize{false};

    auto extracted_mol =
        ExtractMolFragment(mol, selected_atom_indices, sanitize);
    extracted_mol->setProp(HELM_MODEL, true);
    return extracted_mol;
}

} // namespace rdkit_extensions
} // namespace schrodinger
