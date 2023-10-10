#include "schrodinger/rdkit_extensions/coord_utils.h"

#include <Geometry/point.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

unsigned int compute2DCoords(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& frozen_ids)
{

    RDDepict::Compute2DCoordParameters params;
    params.forceRDKit = true;
    params.useRingTemplates = true;

    RDGeom::INT_POINT2D_MAP coord_map;
    if (!frozen_ids.empty()) {
        auto& conf = mol.getConformer();
        for (auto idx : frozen_ids) {
            coord_map[idx] = conf.getAtomPos(idx);
        }
        params.coordMap = &coord_map;
    }

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
