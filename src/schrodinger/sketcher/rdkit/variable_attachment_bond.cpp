#include "schrodinger/sketcher/rdkit/variable_attachment_bond.h"

#include <algorithm>

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/subset.h"

namespace schrodinger
{
namespace sketcher
{

void fix_variable_attachment_bond_coordinates(RDKit::ROMol& mol)
{
    const RDGeom::Point3D X_AXIS(1, 0, 0);
    auto& conf = mol.getConformer();

    for (auto* bond : mol.bonds()) {
        auto variable_attachment_atoms = get_variable_attachment_atoms(bond);
        if (variable_attachment_atoms.empty()) {
            // this isn't a variable attachment bond
            continue;
        }
        RDKit::Atom* dummy_atom = bond->getBeginAtom();
        RDKit::Atom* real_atom = bond->getEndAtom();
        if (real_atom->getAtomicNum() ==
            rdkit_extensions::DUMMY_ATOMIC_NUMBER) {
            std::swap(dummy_atom, real_atom);
        }
        const auto real_atom_pos = conf.getAtomPos(real_atom->getIdx());

        // Calculate new coordinates for the variable attachment bond
        RDGeom::Point3D new_dummy_pos, new_real_atom_pos;
        try {
            std::tie(new_dummy_pos, new_real_atom_pos) =
                get_coordinates_for_variable_attachment_bond(
                    mol, variable_attachment_atoms);
        } catch (variable_attachment_bond_error&) {
            // This variable attachment bond isn't valid (it has either to few
            // variable attachment atoms or the variable attachment atoms aren't
            // involved in any bonds) so we can't fix it
            continue;
        }
        auto new_dummy_angle =
            (new_dummy_pos - new_real_atom_pos).angleTo(X_AXIS);

        // Determine the ideal angle from the real atom to the dummy
        std::vector<RDGeom::Point3D> neighbor_rel_coords;
        for (auto* neighbor : mol.atomNeighbors(real_atom)) {
            if (neighbor == dummy_atom) {
                continue;
            }
            auto& neighbor_pos = conf.getAtomPos(neighbor->getIdx());
            neighbor_rel_coords.push_back(neighbor_pos - real_atom_pos);
        }
        auto ideal_dummy_coordinates =
            best_placing_around_origin(neighbor_rel_coords);
        auto ideal_dummy_angle = ideal_dummy_coordinates.angleTo(X_AXIS);

        auto rotate_by = new_dummy_angle - ideal_dummy_angle;
        auto translate_by = new_real_atom_pos - real_atom_pos;

        // move the variable attachment bond to its new position
        conf.setAtomPos(dummy_atom->getIdx(), new_dummy_pos);
        conf.setAtomPos(real_atom->getIdx(), new_real_atom_pos);

        // Translate any atoms connected to the real atom, and then rotate them
        // about the real atom so that the dummy atom is at the ideal angle.
        auto [connected_atoms, connected_bonds] =
            get_connected_atoms_and_bonds(real_atom);
        for (auto* cur_atom : connected_atoms) {
            if (cur_atom == dummy_atom || cur_atom == real_atom) {
                continue;
            }
            auto& cur_atom_pos = conf.getAtomPos(cur_atom->getIdx());
            cur_atom_pos += translate_by;
            cur_atom_pos = rotate_point_radians(cur_atom_pos, new_real_atom_pos,
                                                rotate_by);
        }
    }
}
} // namespace sketcher
} // namespace schrodinger
