#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"

#include <algorithm>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

RDGeom::Point3D compute_centroid(const std::vector<RDGeom::Point3D>& positions)
{
    // Calculate the centroid by averaging the coordinates
    size_t num_atoms = positions.size();
    RDGeom::Point3D centroid;
    for (const auto& coord : positions) {
        centroid += coord;
    }
    // avoid division by zero
    if (num_atoms > 0) {
        centroid /= static_cast<double>(num_atoms);
    }
    return centroid;
}

unsigned int compute2DCoords(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& frozen_ids)
{
    if (isMonomeric(mol)) {
        return compute_monomer_mol_coords(mol);
    }

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

    // If we are generating a new 2D conformer, we don't want
    // to keep the old bond wedging, as it is relative to the old
    // conformation, which is overwritten by the new one.
    for (auto bond : mol.bonds()) {
        if (bond->hasProp(RDKit::common_properties::_MolFileBondStereo)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondStereo);
            continue;
        }
        if (bond->hasProp(RDKit::common_properties::_MolFileBondCfg)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondCfg);
        }
    }

    return RDDepict::compute2DCoords(mol, params);
}

} // namespace rdkit_extensions
} // namespace schrodinger
