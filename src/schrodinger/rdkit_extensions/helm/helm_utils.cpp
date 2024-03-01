// -------------------------------------------------------------------------
// Implements schrodinger::rdkit_extensions:: helper apis for HELM models
//
//
// Copyright Schrodinger LLC, All Rights Reserved.

#include "schrodinger/rdkit_extensions/helm.h"

#include <optional>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include "schrodinger/container/dynamic_bitset_on_bits_wrapper.hpp"

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
} // namespace rdkit_extensions
} // namespace schrodinger
