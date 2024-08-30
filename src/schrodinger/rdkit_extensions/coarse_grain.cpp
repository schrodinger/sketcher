#include "schrodinger/rdkit_extensions/coarse_grain.h"

#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include <boost/functional/hash.hpp>
#include <fmt/format.h>

#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger::rdkit_extensions
{

ChainType to_chain_type(std::string_view chain_type)
{
    if (chain_type == "PEPTIDE") {
        return ChainType::PEPTIDE;
    } else if (chain_type == "RNA") {
        return ChainType::RNA;
    } else if (chain_type == "DNA") {
        return ChainType::DNA;
    } else if (chain_type == "CHEM") {
        return ChainType::CHEM;
    } else {
        throw std::invalid_argument("Invalid chain type");
    }
}

std::string to_string(ChainType chain_type)
{
    switch (chain_type) {
        case ChainType::PEPTIDE:
            return "PEPTIDE";
        case ChainType::RNA:
            return "RNA";
        case ChainType::DNA:
            return "DNA";
        case ChainType::CHEM:
            return "CHEM";
        default:
            throw std::invalid_argument("Invalid chain type");
    }
}

void add_connection(RDKit::RWMol& cg_mol, size_t monomer1, size_t monomer2,
                    const std::string& linkage)
{
    const auto new_total =
        cg_mol.addBond(monomer1, monomer2, ::RDKit::Bond::BondType::SINGLE);
    cg_mol.getBondWithIdx(new_total - 1)->setProp(LINKAGE, linkage);
}

void add_connection(RDKit::RWMol& cg_mol, size_t monomer1, size_t monomer2,
                    ConnectionType connection_type)
{
    switch (connection_type) {
        case ConnectionType::FORWARD:
            add_connection(cg_mol, monomer1, monomer2, BACKBONE_LINKAGE);
            break;
        case ConnectionType::SIDECHAIN:
            add_connection(cg_mol, monomer1, monomer2, BRANCH_LINKAGE);
            break;
    }
}

std::unique_ptr<Monomer> make_monomer(std::string_view name,
                                      std::string_view chain_id,
                                      int residue_number)
{
    auto a = std::make_unique<::RDKit::Atom>();
    std::string n{name};
    a->setProp(ATOM_LABEL, n);
    a->setProp("Name", n);
    a->setProp("smilesSymbol", n);

    // hack to get some level of canonicalization for CG molecules
    static boost::hash<std::string> hasher;
    a->setIsotope(hasher(n));
    a->setNoImplicit(true);

    auto* residue_info = new ::RDKit::AtomPDBResidueInfo();
    residue_info->setResidueNumber(residue_number);
    residue_info->setResidueName(n);
    residue_info->setChainId(std::string{chain_id});
    a->setMonomerInfo(residue_info);
    return a;
}

size_t add_monomer(RDKit::RWMol& cg_mol, std::string_view name,
                   int residue_number, std::string_view chain_id,
                   MonomerType monomer_type)
{
    auto monomer = make_monomer(name, chain_id, residue_number);
    bool update_label = true;
    bool take_ownership = true;
    auto new_index =
        cg_mol.addAtom(monomer.release(), update_label, take_ownership);
    return new_index;
}

size_t add_monomer(RDKit::RWMol& cg_mol, std::string_view name,
                   MonomerType monomer_type)
{
    if (cg_mol.getNumAtoms() == 0) {
        throw std::invalid_argument(
            "No atoms in molecule to determine chain ID");
    }
    const auto* last_monomer = cg_mol.getAtomWithIdx(cg_mol.getNumAtoms() - 1);
    const auto chain_id = get_polymer_id(last_monomer);
    const auto residue_number = get_residue_number(last_monomer) + 1;
    return add_monomer(cg_mol, name, residue_number, chain_id, monomer_type);
}

void set_residue_number(RDKit::Atom* atom, int residue_number)
{
    auto* residue_info =
        static_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    if (residue_info == nullptr) {
        throw std::runtime_error(
            fmt::format("Atom {} does not have residue info", atom->getIdx()));
    }
    residue_info->setResidueNumber(residue_number);
}

bool is_valid_chain(const RDKit::RWMol& cg_mol, std::string_view polymer_id)
{
    // Check that the residue ordering is valid for this polymer. The residues
    // should be in connectivity order
    const auto chain = get_polymer(cg_mol, polymer_id);
    auto last_residue = chain.atoms[0];
    for (size_t i = 1; i < chain.atoms.size(); ++i) {
        const auto bond =
            cg_mol.getBondBetweenAtoms(last_residue, chain.atoms[i]);
        if (bond == nullptr) {
            return false;
        }

        // Bond direction should be in the same order as residues
        if (get_residue_number(bond->getBeginAtom()) >
                get_residue_number(bond->getEndAtom()) &&
            chain.atoms[i] != bond->getEndAtomIdx()) {
            return false;
        }

        auto linkage = bond->getProp<std::string>(LINKAGE);
        if (linkage == BACKBONE_LINKAGE) {
            last_residue = chain.atoms[i];
        } else if (linkage != BRANCH_LINKAGE) {
            return false;
        }
    }
    return true;
}

RDKit::Atom* find_chain_begin(RDKit::RWMol& cg_mol)
{
    // Find the beginning of the chain by starting at an arbirtary atom
    // and following the backbone backwards until the 'source' (beginning of the
    // chain) is found. If the beginning of the chain is in a cycle, then the
    // last discovered atom will be made the beginning of the chain.
    std::vector<bool> seen(cg_mol.getNumAtoms(), 0);
    auto chain_begin = cg_mol.getAtomWithIdx(0);
    bool updated = true;
    while (updated) {
        updated = false;
        seen[chain_begin->getIdx()] = true;
        for (const auto bond : cg_mol.atomBonds(chain_begin)) {
            auto linkage = bond->getProp<std::string>(LINKAGE);
            if (linkage == "R3-R3") {
                // Don't cross cysteine bridges
                continue;
            } else if (seen[bond->getOtherAtom(chain_begin)->getIdx()]) {
                // this is a loop
                continue;
            }
            // Bonds are always oriented so that the beginAtom->endAtom matches
            // the direction of the chain. If this atom is the 'end' of the
            // bond, it is not the beginning of the chain, so follow the parent
            if (bond->getEndAtom() == chain_begin) {
                auto other = bond->getOtherAtom(chain_begin);
                chain_begin = other;
                updated = true;
                break;
            }
        }
    }
    return chain_begin;
}

void order_residues(RDKit::RWMol& cg_mol)
{
    // Currently assumes that all monomers are in the same chain. We will
    // eventually want to order residues on a per-chain basis.

    // Find the beginning of the chain (must be a backbone monomer)
    auto chain_begin = find_chain_begin(cg_mol);

    // Now re-order the residues beginning at chain_begin
    std::vector<RDKit::Atom*> queue;
    std::vector<bool> visited(cg_mol.getNumAtoms(), 0);
    queue.push_back(chain_begin);
    visited[chain_begin->getIdx()] = true;

    int current_res_idx = 1;
    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        set_residue_number(current, current_res_idx);
        ++current_res_idx;

        // When ordering residues, sidechain monomers should come before
        // backbone monomers, which is more specific then a regular BFS
        // ordering. For example: A.B(X)C should be ordered as A, B, X, C
        // instead of A, B, C, X
        for (const auto bond : cg_mol.atomBonds(current)) {
            if (bond->getEndAtom() == current ||
                visited[bond->getOtherAtom(current)->getIdx()]) {
                continue;
            }
            auto linkage = bond->getProp<std::string>(LINKAGE);
            if (linkage == BRANCH_LINKAGE) {
                auto other = bond->getOtherAtom(current);
                set_residue_number(other, current_res_idx);
                ++current_res_idx;
            } else if (linkage == BACKBONE_LINKAGE) {
                queue.push_back(bond->getOtherAtom(current));
            }
        }
    }
}

void assign_chains(RDKit::RWMol& cg_mol)
{
    cg_mol.setProp("HELM_MODEL", true);

    // Currently, order_residues only works when there is a single chain
    auto chain_ids = get_polymer_ids(cg_mol);
    if (chain_ids.size() == 1 && !is_valid_chain(cg_mol, chain_ids[0])) {
        order_residues(cg_mol);
    }

    // Determine and mark the 'connection bonds'
    if (!cg_mol.getRingInfo()->isInitialized()) {
        ::RDKit::MolOps::findSSSR(cg_mol);
    }
    // get atom rings that belong to a single polymer
    const auto& bnd_rings = cg_mol.getRingInfo()->bondRings();

    for (const auto& ring : bnd_rings) {
        if (std::all_of(ring.begin(), ring.end(), [&](const auto& idx) {
                auto begin_at = cg_mol.getBondWithIdx(idx)->getBeginAtom();
                auto end_at = cg_mol.getBondWithIdx(idx)->getEndAtom();
                return get_polymer_id(begin_at) == get_polymer_id(end_at);
            })) {
            // break this ring -- find bond between residue #s with largest
            // difference
            unsigned int connection_bond = cg_mol.getNumBonds();
            int max_diff = 0;
            for (const auto& idx : ring) {
                const auto bond = cg_mol.getBondWithIdx(idx);
                auto begin_at = bond->getBeginAtom();
                auto end_at = bond->getEndAtom();
                int diff =
                    get_residue_number(begin_at) - get_residue_number(end_at);
                if (std::abs(diff) > max_diff) {
                    connection_bond = bond->getIdx();
                    max_diff = std::abs(diff);
                }
            }
            if (connection_bond != cg_mol.getNumBonds()) {
                std::string linkage = cg_mol.getBondWithIdx(connection_bond)
                                          ->getProp<std::string>(LINKAGE);
                cg_mol.getBondWithIdx(connection_bond)
                    ->setProp(CUSTOM_BOND, linkage);
            } else {
                // temporary error handling
                throw std::runtime_error("Could not find connection bond");
            }
        }
    }
    for (auto bond : cg_mol.bonds()) {
        if (get_polymer_id(bond->getBeginAtom()) !=
            get_polymer_id(bond->getEndAtom())) {
            std::string linkage = bond->getProp<std::string>(LINKAGE);
            bond->setProp(CUSTOM_BOND, linkage);
        }
    }
}

} // namespace schrodinger::rdkit_extensions
