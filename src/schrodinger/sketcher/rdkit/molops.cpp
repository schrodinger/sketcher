#include "schrodinger/sketcher/rdkit/molops.h"

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include <boost/algorithm/string.hpp>

#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/sgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

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

    rdkit_extensions::assign_CIP_labels(*mol);

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
    // update the brackets for the sgroups
    for (auto& sgroup : getSubstanceGroups(mol)) {
        rdkit_extensions::update_s_group_brackets(sgroup);
    }
}

} // namespace sketcher
} // namespace schrodinger
