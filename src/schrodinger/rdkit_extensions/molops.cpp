/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: miscellaneous mol operations
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/molops.h"

#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MolOps.h>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"

namespace schrodinger
{
namespace rdkit_extensions
{

void apply_sanitization(RDKit::RWMol& mol, Sanitization sanitization)
{
    RDLog::LogStateSetter silence_rdkit_logging;

    unsigned failed_op{0};
    int ops = RDKit::MolOps::SANITIZE_ALL;
    if (sanitization == Sanitization::PARTIAL) {
        ops = RDKit::MolOps::SANITIZE_SYMMRINGS |
              RDKit::MolOps::SANITIZE_SETCONJUGATION |
              RDKit::MolOps::SANITIZE_SETHYBRIDIZATION |
              RDKit::MolOps::SANITIZE_ADJUSTHS;
    }
    RDKit::MolOps::sanitizeMol(mol, failed_op, ops);
    // Regardless of sanitization level, ensure property cache is updated
    mol.updatePropertyCache(false);
}

void addHs(RDKit::RWMol& mol, std::vector<unsigned> atom_ids)
{
    // If atom_ids is empty, add Hs to all atoms
    auto only_on_atoms = atom_ids.empty() ? nullptr : &atom_ids;

    auto initial_num_atoms = mol.getNumAtoms();
    bool explicit_only = false;
    bool add_coords = false;
    RDKit::MolOps::addHs(mol, explicit_only, add_coords, only_on_atoms);

    // Workaround for RDKit Issue #7123:
    // addHs sets the NoImplicit flag to true for all atoms; reset it to false
    for (auto idx : atom_ids) {
        mol.getAtomWithIdx(idx)->setNoImplicit(false);
    }

    // Explicitly add 2D coordinates to the new hydrogens; ids are guaranteed to
    // be from the previous number of atoms to the current number of atoms
    std::vector<unsigned int> frozen_ids;
    for (unsigned int idx = 0; idx < initial_num_atoms; ++idx) {
        frozen_ids.push_back(idx);
    }
    rdkit_extensions::compute2DCoords(mol, frozen_ids);
}

void removeHs(RDKit::RWMol& rdk_mol)
{
    RDKit::MolOps::RemoveHsParameters ps;

    // We always remove H on queries; for sketcher import, all atoms are created
    // as QueryAtoms as that they might be changed into queries later on; for
    // conversion from 3D Structure, there is no way to create queries.
    ps.removeWithQuery = true;

    // Disable displaying warnings
    ps.showWarnings = false;

    bool sanitize = false;
    RDKit::MolOps::removeHs(rdk_mol, ps, sanitize);

    rdk_mol.updatePropertyCache(false);
}

void removeHs(RDKit::RWMol& rdk_mol, std::vector<unsigned> atom_ids)
{
    if (atom_ids.empty()) {
        return;
    }

    // Augment atom ids with the ids of the hydrogens attached to them.
    // We don't care about duplicates, because we will sort ant uniquify
    // the list later on.
    const auto num_atom_ids = atom_ids.size();
    for (unsigned i = 0; i < num_atom_ids; ++i) {
        const auto atom = rdk_mol.getAtomWithIdx(atom_ids[i]);
        for (auto nbr : rdk_mol.atomNeighbors(atom)) {
            if (nbr->getAtomicNum() == 1 && nbr->getIsotope() == 0) {
                atom_ids.push_back(nbr->getIdx());
            }
        }
    }

    constexpr int h_protection_mark = 1000;

    // Sort and uniquify the list of atom ids.
    std::sort(atom_ids.begin(), atom_ids.end());
    atom_ids.erase(std::unique(atom_ids.begin(), atom_ids.end()),
                   atom_ids.end());

    auto remove_id_itr = atom_ids.begin();
    std::vector<RDKit::Atom*> protected_atoms;
    const auto input_num_atoms = rdk_mol.getNumAtoms();
    for (unsigned i = 0; i < input_num_atoms; ++i) {
        if (i == *remove_id_itr) {
            ++remove_id_itr;
            continue;
        }
        auto atom = rdk_mol.getAtomWithIdx(i);
        if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
            atom->setIsotope(atom->getIsotope() + h_protection_mark);
            protected_atoms.push_back(atom);
        }
    }

    rdkit_extensions::removeHs(rdk_mol);

    for (auto atom : protected_atoms) {
        atom->setIsotope(atom->getIsotope() - h_protection_mark);
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
