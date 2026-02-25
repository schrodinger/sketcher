#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"
#include <set>

namespace RDKit
{
class ROMol;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{
namespace rdkit_extensions
{
// The RDKit property containing the size of the monomer's graphics item, as
// measured in scene units
const std::string MONOMER_ITEM_SIZE{"monomerItemSize"};

/**
 * Information about a turn in a snaking or coiling chain layout.
 */
struct TurnInfo {
    unsigned int position; // Monomer index where the turn begins (exclusive end
                           // of segment)
    unsigned int
        size; // Number of residues consumed in the turn (0 = instant turn)
    bool downward =
        true; // Whether the turn goes downward (true) or upward (false)
    int y_size_difference =
        0; // this turn will occupy space as if it had y_size_difference more
           // monomers than it actually has. This is used to make sure that
           // turns with different number of monomers still occupy the same
           // vertical space and can be aligned with each other
};

// Structure to hold CUSTOM_BOND information
struct CustomBondInfo {
    unsigned int monomer_i; // First monomer position in chain
    unsigned int monomer_j; // Second monomer position (j > i)
};

/**
 * @param monomer_mol monomeric molecule
 * @return The id of the conformer added to the molecule with the computed
 * coordinates
 */
unsigned int RDKIT_EXTENSIONS_API
compute_monomer_mol_coords(RDKit::ROMol& monomer_mol);

/**
 * Scores the entire coiling layout based on all custom bonds. The lower the
 * score, the better, with 0 being a perfect score where all custom bonds
 * connect residues in corresponding positions at exactly 1 full coil distance.
 *
 * @param polymer The polymer molecule
 * @param turns Vector of turn information defining the layout
 * @param custom_bonds Vector of custom bonds to score
 * @return Overall layout score. Lower is better, with 0 being a perfect score
 */
float RDKIT_EXTENSIONS_API score_coiling_layout(
    const RDKit::ROMol& polymer, const std::vector<TurnInfo>& turns,
    const std::vector<CustomBondInfo>& custom_bonds);

/**
 * @return a float number representing the position of the monomer within the
 * coiling layout.
 * The integer part is the section number (0 for first
 * segment, 1 for first turn, 2 for second segment, etc.), the fractional
 * part is the position within the segment/turn, normalized by the
 * segment/turn size.
 */
float RDKIT_EXTENSIONS_API get_position_in_coils_float(
    unsigned int monomer_idx, const RDKit::ROMol& polymer,
    const std::vector<TurnInfo>& turns);

/**
 * Resize the monomer at the given index to the new size by moving other
 * monomers accordingly
 */
void RDKIT_EXTENSIONS_API
resize_monomers(RDKit::ROMol& monomer_mol,
                std::unordered_map<int, RDGeom::Point3D> monomer_sizes);

/**
 * Check whether any bonds in the given monomer mol cross each other or are
 * closer than a threshold
 */
bool RDKIT_EXTENSIONS_API
has_no_bond_crossings(const RDKit::ROMol& monomer_mol);

/**
 * Check whether any monomers in the given monomer mol are closer than a
 * threshold
 */
bool RDKIT_EXTENSIONS_API has_no_clashes(const RDKit::ROMol& monomer_mol);

/**
 * Check whether the coordinates of the given monomer mol are valid
 * (no clashes, no bond crossings, no stretched bonds)
 */
bool RDKIT_EXTENSIONS_API
coordinates_are_valid(const RDKit::ROMol& monomer_mol);

/**
 * Check whether two line segments defined by points p1, p2 and q1, q2
 * intersect
 */
bool RDKIT_EXTENSIONS_API segments_intersect(const RDGeom::Point3D& p1,
                                             const RDGeom::Point3D& p2,
                                             const RDGeom::Point3D& q1,
                                             const RDGeom::Point3D& q2);

/**
 * Check whether two adjacent bonds (defined by three points) form an angle
 * smaller than a threshold. This is used for pairs of bonds that share an atom,
 * to make sure that the bonds are not collinear, or nearly so.
 */
bool RDKIT_EXTENSIONS_API adjacent_bonds_are_too_close(
    const RDGeom::Point3D& point1, const RDGeom::Point3D& point2,
    const RDGeom::Point3D& point3);

/**
 * Check whether two bonds are too close to each other or intersect. Note that
 * if the two bonds share an atom this will always return true. Use
 * adjacent_bonds_are_too_close() for that case instead
 */
bool RDKIT_EXTENSIONS_API bonds_are_too_close(const RDGeom::Point3D& begin1_pos,
                                              const RDGeom::Point3D& end1_pos,
                                              const RDGeom::Point3D& begin2_pos,
                                              const RDGeom::Point3D& end2_pos);

/* Compute the distance between two line segments defined by points p1, p2 and
 * q1, q2 */
double RDKIT_EXTENSIONS_API compute_distance_between_segments(
    const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
    const RDGeom::Point3D& q1, const RDGeom::Point3D& q2);

/**
 * Initializes RingInfo for `polymer` with all possible rings including any
 * rings that may have been formed with zero order bonds. We do this by
 * converting all the zero order bonds into single order bonds and then looking
 * for rings. (The bond types are restored to their original values by the end
 * of this method)
 */
void RDKIT_EXTENSIONS_API compute_full_ring_info(const RDKit::ROMol& polymer);

// helper function to determine if a ring is layed out as a regular polygon.
// This is done by checking the consistency of the distances from the
// centroid and the uniformity of the angles between consecutive atoms.
bool RDKIT_EXTENSIONS_API is_geometrically_regular_ring_2d(
    const RDKit::ROMol& mol, const std::vector<int>& ring_atoms,
    double radius_tol = 0.1, double angle_tol = 0.25);

/**
 * Find all monomers connected to start_idx that are not part of the given
 * ring. This is used to propagate displacements when resizing rings.
 */
std::set<int> RDKIT_EXTENSIONS_API find_all_connected_monomers_outside_ring(
    const RDKit::ROMol& mol, int start_idx, const std::vector<int>& ring);

} // namespace rdkit_extensions
} // namespace schrodinger
