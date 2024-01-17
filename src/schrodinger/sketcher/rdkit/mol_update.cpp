#include "schrodinger/sketcher/rdkit/mol_update.h"

#include <rdkit/CIPLabeler/CIPLabeler.h>
#include <rdkit/CIPLabeler/TooManyNodesException.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/sgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * @return whether the given mol has CFG/stereo properties present that would
 * have been read from MDL molblock input
 */
static bool has_molblock_cfgs(RDKit::ROMol& mol)
{
    auto bonds = mol.bonds();
    return std::any_of(bonds.begin(), bonds.end(), [](auto& bond) {
        return bond->getBondType() == RDKit::Bond::SINGLE &&
               (bond->hasProp(RDKit::common_properties::_MolFileBondCfg) ||
                bond->hasProp(RDKit::common_properties::_MolFileBondStereo));
    });
}

/**
 * Sketcher-specific assignment of CIP labels, which bypasses the assignment
 * after a certain number of cycles or certain exceptions are hit. For this
 * function to work properly, the input mol must have had stereo assigned.
 */
static void assign_CIP_labels(RDKit::RWMol& mol)
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

void prepare_mol(RDKit::RWMol& mol)
{
    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)
    rdkit_extensions::update_2d_coordinates(mol);

    // SHARED-8774: Deal with chiral flag
    rdkit_extensions::add_enhanced_stereo_to_chiral_atoms(mol);

    // Preserve IDs of enhanced stereo groups from input
    RDKit::forwardStereoGroupIds(mol);

    // Convert parities back to wedges/dashes
    if (has_molblock_cfgs(mol)) {
        RDKit::reapplyMolBlockWedging(mol);
    } else { // Otherwise recalculate chiral bond directions
        rdkit_extensions::wedgeMolBonds(mol, &mol.getConformer());
    }
}

void update_molecule_on_change(RDKit::RWMol& mol)
{
    // Explicitly update the brackets for the sgroups
    rdkit_extensions::update_s_group_brackets(mol);

    // Apply a limited sanitization on molecule to update internal properties
    // while avoiding any changes to the molecule that would alter it to
    // keep the sketcher's input as close to the original as possible.
    rdkit_extensions::apply_sanitization(
        mol, rdkit_extensions::Sanitization::PARTIAL);

    // Convert up/down bonds into parities
    RDKit::MolOps::assignChiralTypesFromBondDirs(mol);
    // Assign stereo from those parities (primarily for bond stereo)
    rdkit_extensions::UseModernStereoPerception use_modern_stereo;
    rdkit_extensions::assign_stereochemistry(mol);
    // Generate R/S/E/Z labels
    assign_CIP_labels(mol);
}

} // namespace sketcher
} // namespace schrodinger
