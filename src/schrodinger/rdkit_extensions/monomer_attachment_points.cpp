// -------------------------------------------------------------------------
// Helpers for attachment points on HELM monomers and their atomistic
// definitions. Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/monomer_attachment_points.h"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

namespace schrodinger::rdkit_extensions
{

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
    const auto* residue = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
        monomer->getMonomerInfo());
    if (residue == nullptr) {
        throw std::invalid_argument("Monomer has no PDB residue info");
    }
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

} // namespace schrodinger::rdkit_extensions
