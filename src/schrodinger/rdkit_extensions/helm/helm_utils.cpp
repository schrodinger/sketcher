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
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/dynamic_bitset_on_bits_wrapper.hpp"
#include "schrodinger/rdkit_extensions/molops.h"

using ::RDKit::RWMol;

namespace schrodinger
{
namespace rdkit_extensions
{

[[nodiscard]] bool isMonomeric(const RDKit::ROMol& mol)
{
    return mol.hasProp(HELM_MODEL);
}

static void check_if_input_is_atomistic(const RDKit::ROMol& mol)
{
    if (!isMonomeric(mol)) {
        throw std::invalid_argument("Atomistic mols are currently unsupported");
    }
}

[[nodiscard]] bool is_dummy_atom(const ::RDKit::Atom* atom)
{
    return atom->hasProp(REPETITION_DUMMY_ID);
}

bool is_dummy_bond(const ::RDKit::Bond* bond)
{
    return is_dummy_atom(bond->getBeginAtom()) ||
           is_dummy_atom(bond->getEndAtom());
}

[[nodiscard]] std::vector<unsigned int>
get_atoms_in_polymer_chain(const RDKit::ROMol& mol, std::string_view polymer_id)
{
    check_if_input_is_atomistic(mol);

    boost::dynamic_bitset<> selected_atoms(mol.getNumAtoms());
    for (const auto atom : mol.atoms()) {
        std::string monomer_polymer_id;
        if (atom->getPropIfPresent(REPETITION_DUMMY_ID, monomer_polymer_id)) {
            if (monomer_polymer_id == polymer_id) {
                selected_atoms.set(atom->getIdx());
            }
        } else if (get_polymer_id(atom) == polymer_id) {
            selected_atoms.set(atom->getIdx());
        }
    }

    schrodinger::dynamic_bitset_on_bits_wrapper on_bits(selected_atoms);
    return {on_bits.begin(), on_bits.end()};
}

[[nodiscard]] std::vector<unsigned int>
get_atoms_in_polymer_chains(const RDKit::ROMol& mol,
                            const std::vector<std::string_view>& polymer_ids)
{
    check_if_input_is_atomistic(mol);

    boost::dynamic_bitset<> selected_atoms(mol.getNumAtoms());
    for (const auto atom : mol.atoms()) {
        std::string polymer_id;
        if (atom->getPropIfPresent(REPETITION_DUMMY_ID, polymer_id)) {
            if (std::find(polymer_ids.begin(), polymer_ids.end(), polymer_id) !=
                polymer_ids.end()) {
                selected_atoms.set(atom->getIdx());
            }
        } else if (std::find(polymer_ids.begin(), polymer_ids.end(),
                             get_polymer_id(atom)) != polymer_ids.end()) {
            selected_atoms.set(atom->getIdx());
        }
    }

    schrodinger::dynamic_bitset_on_bits_wrapper on_bits(selected_atoms);
    return {on_bits.begin(), on_bits.end()};
}

static void
check_for_polymer_groups_and_extended_annotations(const ::RDKit::ROMol& mol)
{
    if (const auto sgroup = get_supplementary_info(mol); sgroup) {
        std::vector<std::string> datafields;
        static const std::string DATAFIELDS{"DATAFIELDS"};
        if (!sgroup->getPropIfPresent(DATAFIELDS, datafields) ||
            datafields.size() < 2u) {
            return;
        }

        if (!datafields[0].empty()) {
            throw std::invalid_argument(
                "Polymer groups are currently unsupported");
        }

        if (!datafields[1].empty()) {
            throw std::invalid_argument(
                "Extended annotations are currently unsupported");
        }
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

[[nodiscard]] std::string get_polymer_id(const ::RDKit::Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getChainId();
}

[[nodiscard]] unsigned int get_residue_number(const ::RDKit::Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getResidueNumber();
}

[[nodiscard]] std::vector<unsigned int>
get_connections(const ::RDKit::ROMol& monomer_mol)
{
    // TODO: check that monomer_mol is a valid monomer mol with chains and
    // custom bonds assigned
    std::vector<unsigned int> connections;

    // any bond with 'customBond'
    for (const auto& bond : monomer_mol.bonds()) {
        if (bond->hasProp(CUSTOM_BOND)) {
            connections.push_back(bond->getIdx());
        }
    }

    std::sort(connections.begin(), connections.end());
    return connections;
}

std::vector<std::string> get_polymer_ids(const RDKit::ROMol& monomer_mol)
{
    std::vector<std::string> polymer_ids;
    for (auto atom : monomer_mol.atoms()) {
        if (is_dummy_atom(atom)) { // query/repeated monomer
            continue;
        }
        auto id = get_polymer_id(atom);
        // in vector to preseve order of polymers
        if (std::find(polymer_ids.begin(), polymer_ids.end(), id) ==
            polymer_ids.end()) {
            polymer_ids.push_back(id);
        }
    }
    return polymer_ids;
}

Chain get_polymer(const RDKit::ROMol& monomer_mol, std::string_view polymer_id)
{
    std::vector<unsigned int> atoms;
    for (auto atom : monomer_mol.atoms()) {
        if (is_dummy_atom(atom)) {
            continue;
        }
        if (get_polymer_id(atom) == polymer_id) {
            atoms.push_back(atom->getIdx());
        }
    }
    // Sort by get_residue_num
    std::sort(atoms.begin(), atoms.end(),
              [&monomer_mol](unsigned int a, unsigned int b) {
                  return get_residue_number(monomer_mol.getAtomWithIdx(a)) <
                         get_residue_number(monomer_mol.getAtomWithIdx(b));
              });
    std::vector<unsigned int> bonds;
    for (auto bond : monomer_mol.bonds()) {
        if (is_dummy_bond(bond)) {
            continue;
        }
        if (get_polymer_id(bond->getBeginAtom()) == polymer_id &&
            get_polymer_id(bond->getEndAtom()) == polymer_id) {
            bonds.push_back(bond->getIdx());
        }
    }

    std::string annotation{};
    for (const auto& sg : ::RDKit::getSubstanceGroups(monomer_mol)) {
        if (!is_polymer_annotation_s_group(sg)) {
            continue;
        }
        if (sg.getProp<std::string>("ID") == polymer_id) {
            annotation = sg.getProp<std::string>(ANNOTATION);
            break;
        }
    }
    return {atoms, bonds, annotation};
}

[[nodiscard]] const RDKit::SubstanceGroup*
get_supplementary_info(const RDKit::ROMol& mol)
{
    for (auto& sgroup : ::RDKit::getSubstanceGroups(mol)) {
        if (is_supplementary_information_s_group(sgroup)) {
            return &sgroup;
        }
    }
    return nullptr;
}

[[nodiscard]] bool
is_supplementary_information_s_group(const RDKit::SubstanceGroup& sgroup)
{
    std::string type;
    std::string fieldname;
    static const std::string TYPE{"TYPE"};
    static const std::string FIELDNAME{"FIELDNAME"};
    if (!(sgroup.getPropIfPresent(TYPE, type) && type == "DAT" &&
          sgroup.getPropIfPresent(FIELDNAME, fieldname) &&
          fieldname == SUPPLEMENTARY_INFORMATION)) {
        return false;
    }

    std::vector<std::string> datafields;
    static const std::string DATAFIELDS{"DATAFIELDS"};
    return sgroup.getPropIfPresent(DATAFIELDS, datafields) &&
           datafields.size() >= 2u;
}

[[nodiscard]] bool
is_polymer_annotation_s_group(const RDKit::SubstanceGroup& sgroup)
{
    std::string type;
    return sgroup.getPropIfPresent("TYPE", type) && type == "COP" &&
           sgroup.hasProp(ANNOTATION) && sgroup.hasProp("ID");
}

} // namespace rdkit_extensions
} // namespace schrodinger
