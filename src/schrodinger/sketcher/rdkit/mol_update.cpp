#include "schrodinger/sketcher/rdkit/mol_update.h"

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>

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

void prepare_mol(RDKit::RWMol& mol)
{
    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)
    rdkit_extensions::update_2d_coordinates(mol);

    // SHARED-8774: Deal with chiral flag
    rdkit_extensions::add_enhanced_stereo_to_chiral_atoms(mol);

    // Preserve IDs of enhanced stereo groups
    RDKit::forwardStereoGroupIds(mol);

    // Convert parities back to wedges/dashes
    if (has_molblock_cfgs(mol)) {
        RDKit::reapplyMolBlockWedging(mol);
    } else { // Otherwise recalculate chiral bond directions
        RDKit::WedgeMolBonds(mol, &mol.getConformer());
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
    rdkit_extensions::assign_CIP_labels(mol);
}

} // namespace sketcher
} // namespace schrodinger
