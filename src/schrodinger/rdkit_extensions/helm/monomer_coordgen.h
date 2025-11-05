#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @param monomer_mol monomeric molecule
 * @return The id of the conformer added to the molecule with the computed
 * coordinates
 */
unsigned int RDKIT_EXTENSIONS_API
compute_monomer_mol_coords(RDKit::ROMol& monomer_mol);

/**
 * Check whether any bonds in the given monomer mol cross each other or are
 * closer than a threshold
 */
bool has_no_bond_crossings(const RDKit::ROMol& monomer_mol);

/**
 * Check whether any monomers in the given monomer mol are closer than a
 * threshold
 */
bool has_no_clashes(const RDKit::ROMol& monomer_mol);

/**
 * Check whether two line segments defined by points p1, p2 and q1, q2
 * intersect
 */
bool segments_intersect(const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
                        const RDGeom::Point3D& q1, const RDGeom::Point3D& q2);

/* Compute the distance between two line segments defined by points p1, p2 and
 * q1, q2 */
double compute_distance_between_segments(const RDGeom::Point3D& p1,
                                         const RDGeom::Point3D& p2,
                                         const RDGeom::Point3D& q1,
                                         const RDGeom::Point3D& q2);

} // namespace rdkit_extensions
} // namespace schrodinger
