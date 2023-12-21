#include "schrodinger/sketcher/rdkit/subset.h"
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/MolOps.h>

namespace schrodinger
{
namespace sketcher
{

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
