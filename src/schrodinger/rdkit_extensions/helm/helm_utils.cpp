// -------------------------------------------------------------------------
// Implements schrodinger::rdkit_extensions:: helper apis for HELM models
//
//
// Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/helm.h"

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cctype>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/dynamic_bitset_on_bits_wrapper.h"
#include "schrodinger/rdkit_extensions/attachment_points.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

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
get_atoms_in_polymer_chain(const RDKit::ROMol& mol,
                           const std::string& polymer_id)
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
                            const std::vector<std::string>& polymer_ids)
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

[[nodiscard]] boost::shared_ptr<RDKit::RWMol>
extract_helm_polymers(const RDKit::ROMol& mol,
                      const std::vector<std::string>& polymer_ids)
{
    check_if_input_is_atomistic(mol);

    auto selected_atom_indices = get_atoms_in_polymer_chains(mol, polymer_ids);
    constexpr bool sanitize{false};

    auto extracted_mol =
        ExtractMolFragment(mol, selected_atom_indices, sanitize);
    extracted_mol->setProp(HELM_MODEL, true);

    // Remove supplementary information S-groups, which contain polymer groups
    // and extended annotations that are not relevant to extracted polymers.
    // Other S-groups (like polymer annotations) are preserved.
    auto& sgroups = RDKit::getSubstanceGroups(*extracted_mol);
    sgroups.erase(std::remove_if(sgroups.begin(), sgroups.end(),
                                 is_supplementary_information_s_group),
                  sgroups.end());
    return extracted_mol;
}

[[nodiscard]] boost::shared_ptr<RDKit::RWMol>
extract_helm_polymer(const RDKit::ROMol& mol, const std::string& polymer_id)
{
    return extract_helm_polymers(mol, {polymer_id});
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

std::optional<unsigned int>
get_attachment_point_num_from_atom_label(const RDKit::Atom& atom)
{
    std::string label;
    atom.getPropIfPresent(RDKit::common_properties::atomLabel, label);
    if (label.starts_with("_R")) {
        try {
            const int number = std::stoi(label.substr(2));
            if (number > 0) {
                return static_cast<unsigned int>(number);
            }
        } catch (const std::logic_error&) {
            // The label is not a numbered attachment point.
        }
    }
    return std::nullopt;
}

std::optional<unsigned int> get_attachment_point_num(const RDKit::Atom& atom)
{
    unsigned int number = 0;
    if (atom.getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                              number) &&
        number > 0) {
        return number;
    }
    if (const auto label_number =
            get_attachment_point_num_from_atom_label(atom)) {
        return label_number;
    }
    if (atom.getAtomicNum() == 0 && atom.getIsotope() > 0) {
        return atom.getIsotope();
    }
    return std::nullopt;
}

std::vector<std::pair<int, std::string>>
get_attachment_points_for_atomistic_mol(const RDKit::ROMol& mol)
{
    std::vector<std::pair<int, std::string>> points;
    for (const auto* atom : mol.atoms()) {
        const auto number = get_attachment_point_num(*atom);
        if (!number) {
            continue;
        }
        const RDKit::Atom* heavy_atom = atom;
        if (atom->getAtomicNum() <= 1) {
            heavy_atom = nullptr;
            for (const auto* neighbor : mol.atomNeighbors(atom)) {
                if (neighbor->getAtomicNum() > 1) {
                    heavy_atom = neighbor;
                    break;
                }
            }
        }
        if (heavy_atom != nullptr) {
            points.emplace_back(*number, heavy_atom->getSymbol());
        }
    }
    std::ranges::sort(points);
    return points;
}

std::vector<std::pair<int, std::string>>
get_attachment_points_for_smiles(const std::string& smiles)
{
    const auto mol = to_rdkit(smiles, Format::EXTENDED_SMILES);
    return get_attachment_points_for_atomistic_mol(*mol);
}

std::vector<std::pair<int, std::string>>
get_attachment_points_for_res(const std::string& resname,
                              const ChainType chain_type)
{
    const auto smiles =
        MonomerDatabase::instance().getMonomerSmiles(resname, chain_type);
    if (!smiles) {
        if (chain_type == ChainType::RNA && resname == "N") {
            return {{1, ""}};
        }
        return {{1, ""}, {2, ""}};
    }
    return get_attachment_points_for_smiles(*smiles);
}

std::vector<std::pair<int, std::string>>
get_numbered_attachment_points_for_monomer(const RDKit::Atom* monomer)
{
    bool is_smiles = false;
    monomer->getPropIfPresent(SMILES_MONOMER, is_smiles);
    if (is_smiles) {
        return get_attachment_points_for_smiles(
            monomer->getProp<std::string>(ATOM_LABEL));
    }
    const auto* residue = static_cast<const RDKit::AtomPDBResidueInfo*>(
        monomer->getMonomerInfo());
    return get_attachment_points_for_res(residue->getResidueName(),
                                         getChainType(*monomer));
}

std::optional<std::string>
get_attachment_point_element_for_monomer(const RDKit::Atom* monomer,
                                         unsigned int attachment_point_num)
{
    std::optional<std::string> found_element;
    for (const auto& [number, element] :
         get_numbered_attachment_points_for_monomer(monomer)) {
        if (number >= 0 &&
            static_cast<unsigned int>(number) == attachment_point_num) {
            found_element = element;
        }
    }
    return found_element;
}

static bool is_na_base(const RDKit::AtomPDBResidueInfo* residue)
{
    const auto& name = residue->getResidueName();
    if (name.empty()) {
        return true;
    }
    const auto last = std::tolower(static_cast<unsigned char>(name.back()));
    return last != 'p' && last != 'r';
}

ConnectorType get_connector_type(const RDKit::Bond* bond,
                                 bool is_secondary_connection)
{
    std::string custom_linkage;
    bond->getPropIfPresent(CUSTOM_BOND, custom_linkage);
    if (custom_linkage == HYDROGEN_BOND_LINKAGE) {
        return ConnectorType::H_BOND;
    }

    if (bond->getBeginAtom()->getMonomerInfo() == nullptr ||
        bond->getEndAtom()->getMonomerInfo() == nullptr) {
        throw std::runtime_error("Atom has no monomer info");
    }
    const auto* begin_residue = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
        bond->getBeginAtom()->getMonomerInfo());
    const auto* end_residue = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
        bond->getEndAtom()->getMonomerInfo());
    if (begin_residue == nullptr || end_residue == nullptr) {
        return ConnectorType::CHEM;
    }

    const auto& begin_id = begin_residue->getChainId();
    const auto& end_id = end_residue->getChainId();
    const bool begin_is_peptide = begin_id.starts_with("PEPTIDE");
    const bool end_is_peptide = end_id.starts_with("PEPTIDE");
    const bool begin_is_rna = begin_id.starts_with("RNA");
    const bool end_is_rna = end_id.starts_with("RNA");
    const bool both_are_peptide_or_rna =
        (begin_is_peptide || begin_is_rna) && (end_is_peptide || end_is_rna);
    if (!both_are_peptide_or_rna) {
        return ConnectorType::CHEM;
    }
    if (!begin_is_peptide && !end_is_peptide) {
        return is_na_base(begin_residue) || is_na_base(end_residue)
                   ? ConnectorType::NA_BACKBONE_TO_BASE
                   : ConnectorType::NA_BACKBONE;
    }

    std::string linkage;
    bond->getPropIfPresent(is_secondary_connection ? CUSTOM_BOND : LINKAGE,
                           linkage);
    if (linkage == PEPTIDE_SIDE_CHAIN_LINKAGE) {
        const auto begin_element = get_attachment_point_element_for_monomer(
            bond->getBeginAtom(), PEPTIDE_R3_ATTACHMENT_POINT);
        const auto end_element = get_attachment_point_element_for_monomer(
            bond->getEndAtom(), PEPTIDE_R3_ATTACHMENT_POINT);
        return begin_element.value_or("") == "S" &&
                       end_element.value_or("") == "S"
                   ? ConnectorType::PEPTIDE_DISULFIDE
                   : ConnectorType::PEPTIDE_SIDE_CHAIN;
    }

    bool begin_is_branch = false;
    bool end_is_branch = false;
    bond->getBeginAtom()->getPropIfPresent(BRANCH_MONOMER, begin_is_branch);
    bond->getEndAtom()->getPropIfPresent(BRANCH_MONOMER, end_is_branch);
    return begin_is_branch || end_is_branch ? ConnectorType::PEPTIDE_BRANCHING
                                            : ConnectorType::PEPTIDE_LINEAR;
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

[[nodiscard]] PolymerConnections
get_interpolymer_connections(const ::RDKit::ROMol& monomer_mol)
{
    PolymerConnections connections;
    for (const auto& polymer_id : get_polymer_ids(monomer_mol)) {
        connections.emplace(polymer_id, std::set<std::string>{});
    }

    for (const auto bond_idx : get_connections(monomer_mol)) {
        const auto* bond = monomer_mol.getBondWithIdx(bond_idx);
        const auto begin_polymer_id = get_polymer_id(bond->getBeginAtom());
        const auto end_polymer_id = get_polymer_id(bond->getEndAtom());
        if (begin_polymer_id == end_polymer_id) {
            continue;
        }
        connections[begin_polymer_id].insert(end_polymer_id);
        connections[end_polymer_id].insert(begin_polymer_id);
    }
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

/**
 * @return the value of the prop_name property from rd_prop, if present.
 * Otherwise, returns an empty string.
 */
static std::string get_string_prop_if_present(const RDKit::RDProps& rd_prop,
                                              const std::string& prop_name)
{
    std::string prop_val;
    rd_prop.getPropIfPresent<std::string>(prop_name, prop_val);
    return prop_val;
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
        if (get_string_prop_if_present(sg, "ID") == polymer_id) {
            annotation = get_string_prop_if_present(sg, ANNOTATION);
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
    return sgroup.getPropIfPresent("TYPE", type) && type == "COP";
}

std::string get_extended_annotations(const ::RDKit::ROMol& mol)
{
    std::string extended_annotations;
    if (mol.getPropIfPresent(EXTENDED_ANNOTATIONS, extended_annotations)) {
        return extended_annotations;
    }

    if (auto sgroup = get_supplementary_info(mol); sgroup) {
        std::vector<std::string> supplementary_info;
        if (sgroup->getPropIfPresent("DATAFIELDS", supplementary_info) &&
            supplementary_info.size() > 1) {
            return supplementary_info[1];
        }
    }

    return {};
}

} // namespace rdkit_extensions
} // namespace schrodinger
