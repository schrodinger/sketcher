#include "schrodinger/sketcher/rdkit/stereochemistry.h"

#include <CIPLabeler/CIPLabeler.h>
#include <CIPLabeler/TooManyNodesException.h>
#include <GraphMol/Chirality.h>

#include "schrodinger/rdkit_extensions/molops.h"

namespace schrodinger
{
namespace sketcher
{

QString get_atom_chirality_label(const RDKit::Atom& atom)
{
    std::string chirality;
    auto stereo_info = RDKit::Chirality::detail::getStereoInfo(&atom);

    // non-CIP ranked atoms have no priority, so potentially chiral
    // atoms should always show an 'undefined' label (SKETCH-1825)
    auto has_non_CIP_neighbor = [](const auto& atom, const auto& stereo_info) {
        auto& mol = atom.getOwningMol();
        for (auto atom_idx : stereo_info.controllingAtoms) {
            auto controlling_atom = mol.getAtomWithIdx(atom_idx);
            if (RDKit::getAtomRLabel(controlling_atom) != 0 ||
                schrodinger::rdkit_extensions::is_attachment_point_dummy(
                    *controlling_atom)) {
                return true;
            }
        }
        return false;
    };

    if (has_non_CIP_neighbor(atom, stereo_info)) {
        return QString(""); // SKETCH-1729: Don't show a defined label
    }

    if (int possible = 0;
        !atom.getPropIfPresent(RDKit::common_properties::_CIPCode, chirality) &&
        atom.getPropIfPresent(RDKit::common_properties::_ChiralityPossible,
                              possible) &&
        possible) {
        return "?"; // possible, but not specified
    }
    return QString(chirality.c_str());
}

QString get_bond_stereo_label(const RDKit::Bond& bond)
{
    std::string label;
    // bond stereo can never be undefined because we always
    // have 2d coords, so double bonds will always be either
    // STEREOANY (no label), E/Z or non-stereo capable (no label).
    bond.getPropIfPresent(RDKit::common_properties::_CIPCode, label);
    return QString::fromStdString(label);
}

void assign_CIP_labels(RDKit::RWMol& mol)
{
    try {
        // This number of calculation cycles takes:
        // ~1s on a Linux Intel(R) Xeon(R) W-2123 CPU @ 3.60GHz.
        // ~0.5s on a 2019 Mac Book Pro with a Intel i7 @ 2.6 GHz.
        // ~1.5s using the WASM sketcher on either of these.
        unsigned max_cycles = 1000000;

        RDKit::CIPLabeler::assignCIPLabels(mol, max_cycles);
    } catch (const RDKit::CIPLabeler::MaxIterationsExceeded&) {
        // CIP label calculation "timed out". Some labels will be omitted.
    } catch (const RDKit::CIPLabeler::TooManyNodesException&) {
        // CIP label calculation graph became too big. It's unlikely we hit
        // this, we'll most probably hit the "time out' before. Still,
        // keep any labels we found so far.
    }
}

} // namespace sketcher
} // namespace schrodinger