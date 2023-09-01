#include "schrodinger/rdkit_extensions/coord_utils.h"

#include <Geometry/point.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

unsigned int compute2DCoords(RDKit::ROMol& mol)
{
    RDDepict::Compute2DCoordParameters params;
    params.useRingTemplates = true;
    return RDDepict::compute2DCoords(mol, params);
}

void update_2d_coordinates(RDKit::ROMol& mol)
{
    auto conf = std::find_if_not(mol.beginConformers(), mol.endConformers(),
                                 std::mem_fn(&RDKit::Conformer::is3D));
    if (conf == mol.endConformers()) {
        compute2DCoords(mol);
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
