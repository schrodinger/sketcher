#include "schrodinger/sketcher/rdkit/subset.h"
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/MolOps.h>

namespace schrodinger
{
namespace sketcher
{

bool is_contiguous_region(std::unordered_set<const RDKit::Atom*> atoms,
                          std::unordered_set<const RDKit::Bond*> bonds)
{
    // create a molecule from the atoms and bonds
    RDKit::RWMol new_mol;
    constexpr bool update_label = true;
    constexpr bool take_ownership = true;
    std::map<unsigned int, unsigned int> atom_map;
    for (auto atom : atoms) {
        RDKit::Atom* new_atom = new RDKit::Atom(atom->getAtomicNum());
        unsigned int new_idx =
            new_mol.addAtom(new_atom, update_label, take_ownership);
        atom_map[atom->getIdx()] = new_idx;
    }
    // Add bonds
    for (auto bond : bonds) {
        unsigned int begin = bond->getBeginAtomIdx();
        unsigned int end = bond->getEndAtomIdx();
        if (atom_map.count(begin) && atom_map.count(end)) {
            new_mol.addBond(atom_map[begin], atom_map[end],
                            bond->getBondType());
        }
    }
    // get the fragments of the molecule
    std::vector<int> frags;
    std::vector<std::vector<int>> frags_mol_atom_mapping;
    RDKit::MolOps::getMolFrags(new_mol, /* sanitizeFrags = */ false, &frags,
                               &frags_mol_atom_mapping,
                               /* copyConformers = */ false);
    return frags_mol_atom_mapping.size() == 1;
}

std::pair<std::unordered_set<const RDKit::Atom*>,
          std::unordered_set<const RDKit::Bond*>>
get_connected_atoms_and_bonds(const RDKit::Atom* const atom)
{
    auto& mol = atom->getOwningMol();
    auto atom_idx = atom->getIdx();
    std::vector<int> frags;
    std::vector<std::vector<int>> frags_mol_atom_mapping;
    // note that RDKit treats zero-order and dative bonds as "real" bonds for
    // connectivity purposes
    RDKit::MolOps::getMolFrags(mol, /* sanitizeFrags = */ false, &frags,
                               &frags_mol_atom_mapping,
                               /* copyConformers = */ false);
    int frag_idx = frags[atom_idx];
    auto& connected_atom_idxs = frags_mol_atom_mapping[frag_idx];
    std::unordered_set<const RDKit::Atom*> connected_atoms;
    std::transform(
        connected_atom_idxs.begin(), connected_atom_idxs.end(),
        std::inserter(connected_atoms, connected_atoms.end()),
        [&mol](int atom_idx) { return mol.getAtomWithIdx(atom_idx); });
    std::unordered_set<const RDKit::Bond*> connected_bonds;
    for (auto* atom : connected_atoms) {
        auto atom_bonds = mol.atomBonds(atom);
        connected_bonds.insert(atom_bonds.begin(), atom_bonds.end());
    }
    return {connected_atoms, connected_bonds};
}

bool in_same_fragment(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto& mol = (*atoms.begin())->getOwningMol();
    std::vector<int> frags;
    std::vector<std::vector<int>> frags_mol_atom_mapping;
    // note that RDKit treats zero-order and dative bonds as "real" bonds for
    // connectivity purposes
    RDKit::MolOps::getMolFrags(mol, /* sanitizeFrags = */ false, &frags,
                               &frags_mol_atom_mapping,
                               /* copyConformers = */ false);
    std::unordered_set<int> atom_frags;
    for (auto atom : atoms) {
        auto frag_idx = frags[atom->getIdx()];
        atom_frags.insert(frag_idx);
    }
    return atom_frags.size() == 1;
}

std::unordered_set<const RDKit::Atom*>
get_smaller_substituent_atoms(const RDKit::ROMol& mol, const RDKit::Bond& bond)
{
    // find the substituents: to do so we remove the bond and find the smallest
    // subgraph that contains either of the bond atoms
    auto new_mol = RDKit::RWMol(mol);
    new_mol.removeBond(bond.getBeginAtomIdx(), bond.getEndAtomIdx());
    std::vector<int> frag_map;
    std::vector<std::vector<int>> frags_mol_atom_mapping;
    RDKit::MolOps::getMolFrags(new_mol, /* sanitizeFrags = */ false, &frag_map,
                               &frags_mol_atom_mapping,
                               /* copyConformers = */ false);
    std::vector<std::vector<int>> substituents;
    // if the molecule contains several disconnected fragments, they will be
    // returned by getMolFrags. We need to keep the ones that contain either of
    // the given bond atoms
    for (auto frag : frags_mol_atom_mapping) {
        if (std::find(frag.begin(), frag.end(), bond.getBeginAtomIdx()) !=
                frag.end() ||
            std::find(frag.begin(), frag.end(), bond.getEndAtomIdx()) !=
                frag.end()) {
            substituents.push_back(frag);
        }
    }
    std::unordered_set<const RDKit::Atom*> atoms;
    if (substituents.size() != 2) {
        // this shouldn't happen, unless the bond was part of a ring. Returning
        // an empty set
        return atoms;
    }
    auto smaller_substituent =
        (substituents[0].size() > substituents[1].size() ? substituents[1]
                                                         : substituents[0]);
    for (auto atom_idx : smaller_substituent) {
        atoms.insert(mol.getAtomWithIdx(atom_idx));
    }
    return atoms;
}

} // namespace sketcher
} // namespace schrodinger
