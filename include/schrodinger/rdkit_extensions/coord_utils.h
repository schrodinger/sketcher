#pragma once

#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace RDKit
{
class Atom;
class Conformer;
class ROMol;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @return the centroid of a set of points
 *
 */
RDKIT_EXTENSIONS_API RDGeom::Point3D
compute_centroid(const std::vector<RDGeom::Point3D>& positions);

/**
 * Generate 2D coordinates for both monomer and atomistic displays. Monomeric
 * molecules use monomeric coordgen. For atomistic molecules with complete PDB
 * residue information, first generate a monomer layout with monomeric
 * coordgen, then lay out the atoms to follow it. Other atomistic molecules
 * use RDKit depiction with ring templates, ignoring SetPreferCoordGen.
 *
 * @param mol RDKit molecule to lay out
 * @param frozen_ids Atom indexes whose coordinates should be preserved
 * @param use_monomer_layout Use the monomer layout as a guide.
 * @return ID of the new 2D conformer
 */
RDKIT_EXTENSIONS_API unsigned int
compute2DCoords(RDKit::ROMol& mol,
                const std::vector<unsigned int>& frozen_ids = {},
                bool use_monomer_layout = true);

} // namespace rdkit_extensions
} // namespace schrodinger
