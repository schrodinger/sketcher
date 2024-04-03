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

static Chain* find_chain_including_monomer(RDKit::ROMol& mol, size_t monomer)
{
    for (auto& sgroup : getSubstanceGroups(mol)) {
        if (sgroup.includesAtom(monomer)) {
            return &sgroup;
        }
    }
    return nullptr;
}

void add_connection(RDKit::RWMol& mol, size_t monomer1, size_t monomer2,
                    const std::string& linkage)
{
    const auto new_total =
        mol.addBond(monomer1, monomer2, ::RDKit::Bond::BondType::SINGLE);
    const auto bond_idx = new_total - 1;
    auto chain = find_chain_including_monomer(mol, monomer1);
    if (chain) {
        chain->addBondWithIdx(bond_idx);
    }
    mol.getBondWithIdx(bond_idx)->setProp(LINKAGE, linkage);
}

void add_connection(RDKit::RWMol& mol, size_t monomer1, size_t monomer2,
                    ConnectionType connection_type)
{
    switch (connection_type) {
        case ConnectionType::FORWARD:
            add_connection(mol, monomer1, monomer2, BACKBONE_LINKAGE);
            break;
        case ConnectionType::SIDECHAIN:
            add_connection(mol, monomer1, monomer2, BRANCH_LINKAGE);
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

size_t add_monomer(Chain& chain, std::string_view name,
                   MonomerType monomer_type)
{
    auto chain_id = chain.getProp<std::string>("ID");
    int residue_number = chain.getAtoms().size() + 1;
    auto& cg_mol = static_cast<RDKit::RWMol&>(chain.getOwningMol());
    auto monomer = make_monomer(name, chain_id, residue_number);
    auto new_index = cg_mol.addAtom(monomer.release());
    chain.addAtomWithIdx(new_index);
    return new_index;
}

Chain& add_chain(RDKit::RWMol& mol, std::string_view chain_name)
{
    Chain chain{&mol, "COP"};
    chain.setProp("ID", std::string{chain_name});
    auto idx = addSubstanceGroup(mol, std::move(chain));

    return getSubstanceGroups(mol)[idx];
}
Chain& add_chain(RDKit::RWMol& mol, ChainType chain_type)
{
    // Add chain
    std::string chain_name;
    switch (chain_type) {
        case ChainType::PEPTIDE:
            chain_name = "PEPTIDE";
            break;
        case ChainType::RNA:
            chain_name = "RNA";
            break;
        case ChainType::DNA:
            chain_name = "DNA";
            break;
        case ChainType::CHEM:
            chain_name = "CHEM";
            break;
    }
    auto& chain =
        add_chain(mol, fmt::format("{}{}", chain_name,
                                   getSubstanceGroups(mol).size() + 1));
    return chain;
}

void assign_chains(RDKit::RWMol& mol)
{
    mol.setProp<bool>(HELM_MODEL, true);

    // TODO: Re-assign chains where necessary

    // The HELM writer only writes bonds in the connection section
    // if those bonds do not belong to a substance group. So, here
    // we find and break rings by removing a bond from a substance group.
    // This would get more complicated with more ring-y CG topologies
    for (auto& sg : getSubstanceGroups(mol)) {
        if (sg.getBonds().size() < sg.getAtoms().size()) {
            // no rings
            continue;
        }
        // This is fragile -- remove bond that is not in HELM chain
        // order, this will be the bond that closes the cycle
        auto remove_bond = sg.getBonds().back();
        for (auto& bond_idx : sg.getBonds()) {
            auto bond = mol.getBondWithIdx(bond_idx);
            if (bond->getBeginAtomIdx() + 1 != bond->getEndAtomIdx()) {
                remove_bond = bond_idx;
                break;
            }
        }
        sg.removeBondWithIdx(remove_bond);
    }
}

} // namespace schrodinger::rdkit_extensions
