#include "schrodinger/sketcher/rdkit/molops.h"

#include <GraphMol/CoordGen.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"

namespace schrodinger
{
namespace sketcher
{
boost::shared_ptr<RDKit::RWMol> text_to_mol(const std::string& text,
                                            Format format)
{
    auto mol = rdkit_extensions::to_rdkit(text, format);

    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)
    get_2d_conformer(*mol);

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

} // namespace sketcher
} // namespace schrodinger
