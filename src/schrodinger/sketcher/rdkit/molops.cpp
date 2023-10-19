#include "schrodinger/sketcher/rdkit/molops.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

#include <boost/algorithm/string.hpp>

#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"

namespace schrodinger
{
namespace sketcher
{
boost::shared_ptr<RDKit::RWMol> text_to_mol(const std::string& text,
                                            const Format format)
{
    auto mol = rdkit_extensions::to_rdkit(text, format);

    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)
    rdkit_extensions::update_2d_coordinates(*mol);

    assign_CIP_labels(*mol);

    // SHARED-8774: Deal with chiral flag
    rdkit_extensions::add_enhanced_stereo_to_chiral_atoms(*mol);

    // Preserve IDs of stereo groups
    RDKit::forwardStereoGroupIds(*mol);

    return mol;
}

void update_molecule_metadata(RDKit::ROMol& mol)
{
    mol.updatePropertyCache(false);
    RDKit::MolOps::fastFindRings(mol);
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

} // namespace sketcher
} // namespace schrodinger
