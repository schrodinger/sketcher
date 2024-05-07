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

void assign_chains(RDKit::RWMol& cg_mol)
{
    // TODO: Re-assign chains where necessary
    cg_mol.setProp("HELM_MODEL", true);

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
