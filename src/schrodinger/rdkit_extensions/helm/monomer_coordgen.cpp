// -------------------------------------------------------------------------
// Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/ChemTransforms/MolFragmenter.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>

#ifdef _MSC_VER
// In Boost 1.81, boost::geometry contains a Windows-only
// "#pragma warning ( pop )" that doesn't have a corresponding push. This
// triggers a compiler warning that gets treated as an error, so we need to
// disable that warning temporarily without using push/pop (otherwise Boost's
// extra pop will pop our push, and then our pop will generate the warning).
// This issue has been fixed in Boost 1.87.
#pragma warning(disable : 4193)
#endif
#include <boost/geometry.hpp>
#ifdef _MSC_VER
#pragma warning(default : 4193)
#endif

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace schrodinger
{
namespace rdkit_extensions
{
// empty space gap between a monomer and the one following it in the chain. This
// will be the length of the visibile bond line connecting the two.
constexpr double SIDE_TO_SIDE_DISTANCE = 0.70;
// when a monomer size is not specified or its value is lower than this, this
// value is used as the minimum size
constexpr double MONOMER_MINIMUM_SIZE = 0.80;

// total distance from the center of one monomer to the center of the following
// in the chain
constexpr double MONOMER_BOND_LENGTH =
    SIDE_TO_SIDE_DISTANCE + MONOMER_MINIMUM_SIZE;

// maximum allowed bond length as a multiple of the ideal bond length. Any bond
// longer than this will be considered "stretched"
constexpr double MAX_BOND_STRETCH = 3.0;
constexpr double NON_BACKBONE_BOND_FLEXIBILITY = 2.0;

constexpr double DIST_BETWEEN_MULTIPLE_POLYMERS = 5;
constexpr unsigned int MAX_MONOMERS_PER_SNAKE = 10;
const double PI = boost::math::constants::pi<double>();

constexpr double MONOMER_CLASH_DISTANCE = MONOMER_BOND_LENGTH * 0.25;
constexpr double BOND_CLASH_DISTANCE = MONOMER_CLASH_DISTANCE;
constexpr double MIN_ANGLE_BETWEEN_ADJACENT_BONDS =
    boost::math::constants::pi<double>() / 18; // 10 degrees

constexpr double EPSILON = 1e-6;

// Props used in the fragmented polymer mols to store monomer graph info from
// the monomersitc mol
const std::string BOND_TO{"bondTo"};
const std::string ORIGINAL_INDEX{"originalIndex"};

// Set on the polymer mol
const std::string POLYMER_ID{"polymerID"};

// Replacement atom label for monomers with SMILES strings as their atom label
const std::string SMILES_MONOMER_LABEL{"CX"};

// the direction in which to lay out a monomer chain (from left to right or
// right to left)
enum class ChainDirection { LTR, RTL };

// the general direction in which to lay out branch monomers (above or below the
// chain)
enum class BranchDirection { UP, DOWN };

// the direction in which to place branch monomers relative to their parent
// monomer. This might be differtent from the general branch direction if there
// are clashes.
enum class Direction { N, S, E, W, NW, NE, SW, SE };
enum class PolygonStartSide { LEFT, RIGHT };
static const std::map<Direction, RDGeom::Point3D> DIRECTION_TO_POINT_MAP = {
    {Direction::N, RDGeom::Point3D(0, 1, 0)},
    {Direction::S, RDGeom::Point3D(0, -1, 0)},
    {Direction::E, RDGeom::Point3D(1, 0, 0)},
    {Direction::W, RDGeom::Point3D(-1, 0, 0)},
    {Direction::NE, RDGeom::Point3D(1, 1, 0)},
    {Direction::NW, RDGeom::Point3D(-1, 1, 0)},
    {Direction::SE, RDGeom::Point3D(1, -1, 0)},
    {Direction::SW, RDGeom::Point3D(-1, -1, 0)}};

static RDGeom::Point3D direction_to_point(Direction dir)
{
    auto it = DIRECTION_TO_POINT_MAP.find(dir);
    if (it != DIRECTION_TO_POINT_MAP.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid direction");
}

static void place_monomer_at(RDKit::Conformer& conformer,
                             const RDKit::Atom* monomer_to_place,
                             const RDGeom::Point3D& position,
                             std::unordered_set<int>& placed_monomers_idcs)
{
    auto monomer_idx = monomer_to_place->getIdx();
    conformer.setAtomPos(monomer_idx, position);
    placed_monomers_idcs.insert(monomer_to_place->getProp<int>(ORIGINAL_INDEX));
}

/**
 * Compute the centroid of a set of atom indices in a molecule.
 */
RDGeom::Point3D compute_centroid(const RDKit::ROMol& mol,
                                 const std::vector<int>& atom_indices)
{
    if (atom_indices.empty()) {
        return RDGeom::Point3D(0, 0, 0);
    }
    const auto& conformer = mol.getConformer();
    std::vector<RDGeom::Point3D> positions(atom_indices.size());
    std::transform(atom_indices.begin(), atom_indices.end(), positions.begin(),
                   [&conformer](int idx) { return conformer.getAtomPos(idx); });

    return compute_centroid(positions);
}

/**
 * Check if a bond contains a backbone linkage
 * @return true if the bond is a backbone connection, or if it has more than one
 * connection (i.e., the LINKAGE property contains more than just the
 * CUSTOM_BOND linkage), in which case one is assumed to be a backbone
 * connection.
 */
static bool has_backbone_linkage(const RDKit::Bond* bond)
{
    if (!bond->hasProp(CUSTOM_BOND)) {
        // No custom bond property, this is a regular backbone bond
        return true;
    }

    std::string custom_bond;
    std::string linkage;

    if (!bond->getPropIfPresent(CUSTOM_BOND, custom_bond) ||
        !bond->getPropIfPresent(LINKAGE, linkage)) {
        return false;
    }
    // If LINKAGE and CUSTOM_BOND are different, it means there are additional
    // linkages beyond the custom one: assume that one is a backbone connection
    return linkage != custom_bond;
}

namespace bg = boost::geometry;
using point_t = bg::model::d2::point_xy<double>;
using segment_t = bg::model::segment<point_t>;

static auto default_stop_condition = [](const RDKit::Atom* monomer) {
    return false;
};

/**
 * Get the available directions in which branches for a given monomer can be
 * placed
 * @param bonded_to_parent_polymer true if the monomer has a bond to an already
 * placed polymer
 * @param bonded_to_child_polymer true if the monomer has a bond to an unplaced
 * polymer
 * @param branch_direction the general direction in which branches should be
 * placed
 * @param chain_dir the direction in which the monomer chain is being laid out
 * @param is_last_of_chain true if the monomer is the last in the chain
 * @param is_first_of_chain true if the monomer is the first in the chain
 */
static std::vector<Direction> get_available_directions(
    const bool bonded_to_parent_polymer, const bool bonded_to_child_polymer,
    const BranchDirection branch_direction, const ChainDirection chain_dir,
    const bool is_last_of_chain, const bool is_first_of_chain)
{
    std::vector<Direction> dirs;
    // first try the main branch direction
    if (!bonded_to_child_polymer && branch_direction == BranchDirection::DOWN) {
        dirs.push_back(Direction::S);
    } else if (!bonded_to_parent_polymer &&
               branch_direction == BranchDirection::UP) {
        dirs.push_back(Direction::N);
    }
    // if this is the first or last polymer, we can branch along
    // the chain
    if (is_last_of_chain) {
        dirs.push_back(chain_dir == ChainDirection::LTR ? Direction::E
                                                        : Direction::W);
    }
    if (is_first_of_chain) {
        dirs.push_back(chain_dir == ChainDirection::LTR ? Direction::W
                                                        : Direction::E);
    }
    // as a last resort, try the opposite of the main branch direction
    if (!bonded_to_parent_polymer &&
        branch_direction == BranchDirection::DOWN) {
        dirs.push_back(Direction::N);
    } else if (!bonded_to_child_polymer &&
               branch_direction == BranchDirection::UP) {
        dirs.push_back(Direction::S);
    }
    // add diagonal directions. These are not normally used, but can if
    // necessary
    dirs.push_back(Direction::NE);
    dirs.push_back(Direction::SW);
    dirs.push_back(Direction::NW);
    dirs.push_back(Direction::SE);
    return dirs;
}
/**
 * Lays out a monomer chain in `polymer` by traversing through the monomer graph
 * until either the `stop_condition` is met for all visited monomers or we run
 * out of unplaced monomers reachable from start_monomer. Expects all monomers
 * to have their indices in `placed_monomers_idcs` to indicate if they've
 * already been placed (polymers returned from `break_into_polymers` below
 * already have this prop set on all monomers).
 */
template <typename _Cond = decltype(default_stop_condition)> static void
lay_out_chain(RDKit::ROMol& polymer, const RDKit::Atom* start_monomer,
              std::unordered_set<int>& placed_monomers_idcs,
              const RDGeom::Point3D& start_pos = RDGeom::Point3D(0, 0, 0),
              ChainDirection chain_dir = ChainDirection::LTR,
              BranchDirection branch_direction = BranchDirection::UP,
              _Cond stop_condition = default_stop_condition)
{
    auto& conformer = polymer.getConformer();
    auto x_pos = start_pos.x;
    auto monomer_to_place = start_monomer;
    const RDKit::Atom* last_placed_monomer = nullptr;
    while (monomer_to_place != nullptr) {
        if (stop_condition(monomer_to_place)) {
            break;
        }
        place_monomer_at(conformer, monomer_to_place,
                         RDGeom::Point3D(x_pos, start_pos.y, start_pos.z),
                         placed_monomers_idcs);
        auto monomer_neighbors = polymer.atomNeighbors(monomer_to_place);
        bool bonded_to_parent_polymer = false;
        bool bonded_to_child_polymer = false;
        std::vector<int> neighbor_monomer_idcs;
        monomer_to_place->getPropIfPresent<std::vector<int>>(
            BOND_TO, neighbor_monomer_idcs);
        for (auto neighbor_idx : neighbor_monomer_idcs) {
            if (placed_monomers_idcs.contains(neighbor_idx)) {
                bonded_to_parent_polymer = true;
            } else {
                bonded_to_child_polymer = true;
            }
        }
        auto monomer_to_place_idx = monomer_to_place->getIdx();
        monomer_to_place = nullptr;
        const RDKit::Atom* next_monomer_to_place = nullptr;

        std::vector<const RDKit::Atom*> branches;
        for (auto neighbor : monomer_neighbors) {
            if (placed_monomers_idcs.contains(
                    neighbor->getProp<int>(ORIGINAL_INDEX)) ||
                // Skip if the bond to this neighbor is not a backbone
                // connection)
                !has_backbone_linkage(polymer.getBondBetweenAtoms(
                    neighbor->getIdx(), monomer_to_place_idx))) {
                continue;
            }

            if (!neighbor->getProp<bool>(BRANCH_MONOMER)) {
                next_monomer_to_place = neighbor;
            } else {
                branches.push_back(neighbor);
            }
        }

        std::vector<Direction> available_directions = get_available_directions(
            bonded_to_parent_polymer, bonded_to_child_polymer, branch_direction,
            chain_dir, next_monomer_to_place == nullptr,
            last_placed_monomer == nullptr);
        auto next_available_direction = available_directions.begin();
        if (available_directions.size() < branches.size()) {
            throw std::runtime_error(
                "Not enough available directions to place all branch monomers");
        }

        for (auto branch_monomer : branches) {
            RDGeom::Point3D pos(x_pos, start_pos.y, start_pos.z);
            pos += direction_to_point(*next_available_direction++) *
                   MONOMER_BOND_LENGTH;
            place_monomer_at(conformer, branch_monomer, pos,
                             placed_monomers_idcs);
        }
        last_placed_monomer = monomer_to_place;
        monomer_to_place = next_monomer_to_place;

        x_pos += chain_dir == ChainDirection::LTR ? MONOMER_BOND_LENGTH
                                                  : -MONOMER_BOND_LENGTH;
    }
}

static bool is_nucleic_acid(const RDKit::ROMol& polymer)
{
    auto polymer_id = polymer.getProp<std::string>(POLYMER_ID);
    return boost::starts_with(polymer_id, "DNA") ||
           boost::starts_with(polymer_id, "RNA");
}

void compute_full_ring_info(const RDKit::ROMol& polymer)
{
    std::vector<unsigned int> zob_idxs{};
    for (auto bond : polymer.bonds()) {
        if (bond->getBondType() == RDKit::Bond::BondType::ZERO) {
            zob_idxs.push_back(bond->getIdx());
            bond->setBondType(RDKit::Bond::BondType::SINGLE);
        }
    }

    if (polymer.getRingInfo()->isInitialized()) {
        polymer.getRingInfo()->reset();
    }
    constexpr bool include_dative_bonds = true;
    RDKit::MolOps::findSSSR(polymer, /*res=*/nullptr, include_dative_bonds);

    for (auto bond_idx : zob_idxs) {
        const_cast<RDKit::Bond*>(polymer.getBondWithIdx(bond_idx))
            ->setBondType(RDKit::Bond::BondType::ZERO);
    }
}

/**
 * Generates co-ordinates for a regular n-sided polygon with sides of length
 * `MONOMER_BOND_LENGTH`. The polygon will be oriented so the first edge will be
 * parallel to the y-axis.
 * @param n
 * @param start_side The side of the polygon that the first edge (which is
 * parallel to the y-axis) should be on
 */
static std::vector<RDGeom::Point3D>
get_coords_for_ngon(size_t n,
                    PolygonStartSide start_side = PolygonStartSide::RIGHT)
{
    // radius for n-sided polygon with equal sides
    double radius = MONOMER_BOND_LENGTH / (2 * std::sin(PI / n));
    // rotate to keep first edge parallel with y-axis
    double rotation = -1 * PI / n;

    int x_dir = start_side == PolygonStartSide::LEFT ? -1 : 1;
    std::vector<RDGeom::Point3D> points{};
    for (size_t i = 0; i < n; i++) {
        points.emplace_back(radius * std::cos(2 * PI * i / n + rotation) *
                                x_dir,
                            radius * std::sin(2 * PI * i / n + rotation), 0);
    }
    return points;
}

static RDGeom::Point3D rotate_point(const RDGeom::Point3D& point,
                                    float sin_angle, float cos_angle,
                                    const RDGeom::Point3D& center)
{
    RDGeom::Point3D rotated_pos = point;
    // Translate to origin (rotation center)
    rotated_pos -= center;
    // Rotate
    double x_new = rotated_pos.x * cos_angle - rotated_pos.y * sin_angle;
    double y_new = rotated_pos.x * sin_angle + rotated_pos.y * cos_angle;
    rotated_pos.x = x_new;
    rotated_pos.y = y_new;
    // Translate back
    rotated_pos += center;
    return rotated_pos;
}

/**
 * Rotate all atoms in the conformer by the given angle around the given center.
 * @param conformer The conformer to rotate
 * @param angle Rotation angle in radians (positive = counterclockwise)
 * @param center The center point for rotation
 */
static void rotate_conformer(RDKit::Conformer& conformer, double angle,
                             const RDGeom::Point3D& center)
{
    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);
    for (size_t i = 0; i < conformer.getNumAtoms(); ++i) {
        auto pos = conformer.getAtomPos(i);
        pos = rotate_point(pos, sin_angle, cos_angle, center);
        conformer.setAtomPos(i, pos);
    }
}

/**
 * Orient the ring system based on attachment points to non-ring atoms. This is
 * important for ensuring that the chains come out of the rings in a consistent
 * direction and extend horizontally as much as possible. If there is one
 * attachment point, orient it to the right of the ring's centroid. If there's
 * two attachment points, align them vertically with the rings to the left. If
 * there are more than two attachment points, use the first two (this is
 * somewhat arbitrary but cases with more than two attachment points are
 * expected to be rare and it's not clear what a better strategy would be).
 */
static void orient_ring_system(RDKit::ROMol& polymer)
{
    auto& conformer = polymer.getConformer();

    // Find ring atoms that have bonds to non-ring atoms and store them along
    // with the centroids of their rings
    std::vector<std::pair<unsigned int, RDGeom::Point3D>>
        attachment_atoms_and_centroids;

    for (const auto& ring : polymer.getRingInfo()->atomRings()) {
        // compute centroid of the ring
        auto centroid = compute_centroid(polymer, ring);
        for (auto atom_idx : ring) {
            auto neighbors =
                polymer.atomNeighbors(polymer.getAtomWithIdx(atom_idx));
            for (auto neighbor : neighbors) {
                if (polymer.getRingInfo()->numAtomRings(neighbor->getIdx()) ==
                    0) {
                    attachment_atoms_and_centroids.push_back(
                        {atom_idx, centroid});
                    break;
                }
            }
        }
    }
    // If no attachment atoms found, nothing to orient
    if (attachment_atoms_and_centroids.empty()) {
        return;
    }

    double rotation_angle = 0.0;
    RDGeom::Point3D rotation_center(0, 0, 0);

    if (attachment_atoms_and_centroids.size() == 1) {
        // Single attachment point: rotate so it's to the right of its ring
        // centroid
        auto& [atom_idx, centroid] = attachment_atoms_and_centroids[0];
        auto atom_pos = conformer.getAtomPos(atom_idx);

        // Calculate angle from centroid to attachment point
        double dx = atom_pos.x - centroid.x;
        double dy = atom_pos.y - centroid.y;
        double current_angle = std::atan2(dy, dx);

        // Target angle is 0 (pointing right)
        rotation_angle = -current_angle;
        rotation_center = centroid;

    } else { // at least two attachment points: orient the first two attachment
             // points vertically
        auto& [atom_idx1, centroid1] = attachment_atoms_and_centroids[0];
        auto& [atom_idx2, centroid2] = attachment_atoms_and_centroids[1];
        auto atom1_pos = conformer.getAtomPos(atom_idx1);
        auto atom2_pos = conformer.getAtomPos(atom_idx2);

        // Calculate angle of line between the two attachment points
        double current_angle =
            std::atan2(atom2_pos.y - atom1_pos.y, atom2_pos.x - atom1_pos.x);
        double target_angle = PI / 2; // vertical
        rotation_angle = target_angle - current_angle;

        // Center on midpoint of the two attachment atoms
        rotation_center = (atom1_pos + atom2_pos) * 0.5;

        // If the midpoint of the two centroids would fall to the right, add
        // 180° to the rotation, since we want the rings to the left of the
        // attachment points instead of the right
        auto centroids_midpoint = (centroid1 + centroid2) * 0.5;
        auto rotated_centroid_midpoint =
            rotate_point(centroids_midpoint, std::sin(rotation_angle),
                         std::cos(rotation_angle), rotation_center);
        if (rotated_centroid_midpoint.x > rotation_center.x + EPSILON) {
            rotation_angle += PI;
        }
    }

    // Rotate the whole ring system using the helper function
    rotate_conformer(conformer, rotation_angle, rotation_center);
}

/**
 * Generate coordinates for the ring-systems in `polymer`. The coordinates are
 * generated using RDKit's internal coordinate generation. This is used for
 * cyclic polymers and also as a fallback when the snaking layout fails (e.g.
 * due to incompatible turn constraints from multiple CUSTOM_BONDS). After
 * generating the coordinates, the rings are oriented so that the first
 * CUSTOM_BOND, if any are present, in the rings is vertical. Coordinates will
 * be assigned to all monomers, but only the ones in rings and their immediate
 * neighbors will be marked as placed in `placed_monomers_idcs`, so that the
 * other coordinates can be overwritten later (e.g. to draw straight chains)
 */
static void
generate_coordinates_for_cycles(RDKit::ROMol& polymer,
                                std::unordered_set<int>& placed_monomers_idcs)
{
    // generate ring coordinates using RDKit's small-molecule built-in
    // coordinate generation

    RDDepict::Compute2DCoordParameters params;
    params.forceRDKit = true;
    params.useRingTemplates = false;
    RDDepict::compute2DCoords(polymer, params);
    orient_ring_system(polymer);

    // finalize coordinates of rings and of their immediate neighbors.

    std::set<int> monomers_to_place;
    for (auto ring_atoms : polymer.getRingInfo()->atomRings()) {
        for (auto monomer_idx : ring_atoms) {
            monomers_to_place.insert(monomer_idx);
            auto monomer = polymer.getAtomWithIdx(monomer_idx);
            for (auto neighbor : polymer.atomNeighbors(monomer)) {
                monomers_to_place.insert(neighbor->getIdx());
            }
        }
    }
    // Coordinates need to be scaled from RDKit's default bond length to
    // MONOMER_BOND_LENGTH
    auto scale_factor = MONOMER_BOND_LENGTH / RDDepict::BOND_LEN;
    for (auto monomer_idx : monomers_to_place) {
        auto monomer = polymer.getAtomWithIdx(monomer_idx);
        place_monomer_at(polymer.getConformer(), monomer,
                         polymer.getConformer().getAtomPos(monomer_idx) *
                             scale_factor,
                         placed_monomers_idcs);
    }
}

// Structure to hold turn constraint: turn must be in range (min_pos, max_pos)
struct TurnConstraint {
    unsigned int min_pos; // Turn must be after this position
    unsigned int max_pos; // Turn must be before this position
};

/**
 * Extracts the CUSTOM_BONDs from `polymer` and returns them as a list of
 * CustomBondInfo, sorted by starting position.
 */
static std::vector<CustomBondInfo>
extract_custom_bonds(const RDKit::ROMol& polymer)
{
    std::vector<CustomBondInfo> custom_bonds;
    for (const auto& bond : polymer.bonds()) {
        if (bond->hasProp(CUSTOM_BOND)) {
            auto begin_idx = bond->getBeginAtom()->getIdx();
            auto end_idx = bond->getEndAtom()->getIdx();

            unsigned int i = std::min(begin_idx, end_idx);
            unsigned int j = std::max(begin_idx, end_idx);
            custom_bonds.push_back({i, j});
        }
    }

    // Sort bonds by starting position
    std::sort(
        custom_bonds.begin(), custom_bonds.end(),
        [](const auto& a, const auto& b) { return a.monomer_i < b.monomer_i; });
    return custom_bonds;
}

/**
 * Extracts turn constraints from the CUSTOM_BONDs in `polymer`. Each
 * CUSTOM_BOND creates a constraint that there must be exactly one turn between
 * the two monomers it connects, so they end up in consecutive rows in the
 * snaking layout and can be connected by a single bond. If there are multiple
 * bonds that create impossible layouts (e.g. the two bonds would cross), an
 * empty vector is returned to indicate that the snaking layout is not possible
 * and a fallback layout should be used instead.
 */
std::vector<TurnConstraint>
compute_turn_positions_for_chain(const RDKit::ROMol& polymer)
{
    // Extract all CUSTOM_BONDS
    std::vector<CustomBondInfo> custom_bonds = extract_custom_bonds(polymer);

    // No CUSTOM_BONDS: this should not happens since linear polymers should be
    // handled by a different code path. Return empty constraints so the calling
    // code can use a different strategy
    if (custom_bonds.empty()) {
        return {};
    }

    // Each CUSTOM_BOND requires exactly one turn between i and j
    std::vector<TurnConstraint> constraints;
    for (const auto& bond : custom_bonds) {
        constraints.push_back({bond.monomer_i, bond.monomer_j});
    }

    // If the regions "pinched" by two bonds are one inside the other (i1 < i2 <
    // j2 < j1), they  they share a turn between i2 and j2
    for (size_t idx1 = 0; idx1 < custom_bonds.size(); ++idx1) {
        for (size_t idx2 = idx1 + 1; idx2 < custom_bonds.size(); ++idx2) {
            auto& bond1 = custom_bonds[idx1];
            auto& bond2 = custom_bonds[idx2];

            // check if bond2 is  contained within bond1: i1 < i2 < j2 < j1
            if (bond1.monomer_i < bond2.monomer_i &&
                bond2.monomer_i < bond2.monomer_j &&
                bond2.monomer_j < bond1.monomer_j) {

                //  bond2 is completely contained within bond1 - they will
                //  connect the same pair of consecutive rows The turn must be
                //  between i2 and j1
                constraints.push_back({bond2.monomer_i, bond2.monomer_j});
            } else {
                // bond regions overlap. This means that the two bonds would
                // cross, and the layout is not possible
                return {};
            }
        }
    }

    // Check if all constraints are satisfiable
    // Merge overlapping constraints to find actual turn positions
    std::vector<TurnConstraint> merged;
    for (const auto& constraint : constraints) {
        bool merged_into_existing = false;

        for (auto& existing : merged) {
            // Check if constraints overlap
            if (!(constraint.max_pos <= existing.min_pos ||
                  constraint.min_pos >= existing.max_pos)) {
                // Merge: take intersection of valid ranges
                existing.min_pos =
                    std::max(existing.min_pos, constraint.min_pos);
                existing.max_pos =
                    std::min(existing.max_pos, constraint.max_pos);
                merged_into_existing = true;
                break;
            }
        }

        if (!merged_into_existing) {
            merged.push_back(constraint);
        }
    }

    // Check if any merged constraint is unsatisfiable
    for (const auto& constraint : merged) {
        if (constraint.min_pos >= constraint.max_pos) {
            return {};
        }
    }

    return merged;
}

/**
 * Lays out turn monomers following a polygon arc.
 *
 * @param polymer The polymer molecule
 * @param conformer The conformer to modify
 * @param placed_monomers_idcs Set to track which monomers have been placed
 * @param turn_start_pos Starting position for the turn (end of previous
 * segment)
 * @param segment_end Index where the previous segment ended
 * @param turn_size Number of turn monomers to place
 * @param total_monomers Total number of monomers in the polymer
 * @param chain_dir Direction of the chain segment that was just placed
 * @param downward Whether the turn should go downward (true) or upward (false)
 * @param y_size_difference this turn will occupy space as if it had
 * y_size_difference more monomers than it actually has. This is used to make
 * sure that turns with different number of monomers still occupy the same
 * vertical space and can be aligned with each other
 * @return The position where the next segment should start
 */
static RDGeom::Point3D
lay_out_turn(RDKit::ROMol& polymer, RDKit::Conformer& conformer,
             std::unordered_set<int>& placed_monomers_idcs,
             const RDGeom::Point3D& turn_start_pos, unsigned int segment_end,
             unsigned int turn_size, unsigned int total_monomers,
             ChainDirection chain_dir, bool downward = true,
             int y_size_difference = 0)
{
    if (turn_size == 0) {
        // No turn monomers, just return position one bond length below
        return RDGeom::Point3D(turn_start_pos.x,
                               turn_start_pos.y - MONOMER_BOND_LENGTH, 0.0);
    }

    RDGeom::Point3D current_pos = turn_start_pos;

    // Calculate external angle for regular polygon with (turn_size + 2) * 2
    // vertices
    unsigned int n = 2 * (turn_size + 2);

    float bond_scale = 1.0;
    if (y_size_difference != 0) {
        // scale y coordinates by a factor to make sure that the turn occupies
        // the same vertical space as if it had y_size_difference more monomers.
        // This is the ratio between the apothem of the regular polygon with
        // turn_size monomers and the apothem of the regular polygon with
        // turn_size + y_size_difference monomers
        bond_scale = (2 * std::tan(PI / n)) /
                     (2 * std::tan(PI / (n + 2 * y_size_difference)));
    }

    double external_angle = -2.0 * M_PI / n;

    // For RTL, rotate in opposite direction (clockwise instead of
    // counterclockwise)
    if (chain_dir == ChainDirection::RTL) {
        external_angle = -external_angle;
    }
    // if the turn goes upward instead of downward, also rotate in opposite
    // direction
    if (!downward) {
        external_angle = -external_angle;
    }

    // Start with horizontal direction based on chain direction
    // LTR = 0° (east), RTL = 180° (west)
    double current_angle = (chain_dir == ChainDirection::LTR) ? 0.0 : M_PI;

    // Place each turn monomer by rotating the direction
    for (unsigned int i = 0; i < turn_size; ++i) {
        unsigned int monomer_idx = segment_end + i;
        if (monomer_idx >= total_monomers) {
            break;
        }

        // Rotate by external angle
        current_angle += external_angle;

        // Calculate step using 2D rotation
        double step_x =
            MONOMER_BOND_LENGTH * std::cos(current_angle) * bond_scale;
        double step_y =
            MONOMER_BOND_LENGTH * std::sin(current_angle) * bond_scale;

        current_pos.x += step_x;
        current_pos.y += step_y;

        conformer.setAtomPos(monomer_idx, current_pos);
        placed_monomers_idcs.insert(monomer_idx);
    }

    // Calculate position for next segment: one more rotation and step
    current_angle += external_angle;
    double final_step_x =
        MONOMER_BOND_LENGTH * std::cos(current_angle) * bond_scale;
    double final_step_y =
        MONOMER_BOND_LENGTH * std::sin(current_angle) * bond_scale;

    return RDGeom::Point3D(current_pos.x + final_step_x,
                           current_pos.y + final_step_y, 0.0);
}

/**
 * Computes the turn positions and sizes for a coiling layout based on the the
 * following parameters:
 * @param polymer The polymer molecule
 * @param increment The amount by which the turn size should increase for each
 * subsequent turn. This controls how quickly the coil expands and can be used
 * to approximate an elliptical shape. If it has a negative value the coiling
 * will be reversed, starting with the longest turn and decreasing the turn size
 * for each subsequent turn.
 * @param turn_size The size of the first turn (number of monomers in the turn)
 * @param chain_size The size of the chain segments between turns (number of
 * monomers in each segment). This remains constant for every coil
 * @param first_chain_size The size of the first chain segment (number of
 * monomers in the first segment). This is used to allow the first turn to be
 * placed later or earlier than the regular chain_size, which can be useful to
 * better align the coil with certain CUSTOM_BOND positions.
 * @return A vector of TurnInfo with the position and size and direction
 * (downward or upward) of each turn. If the parameters are invalid an empty
 * vector is returned to indicate failure.
 */
static std::vector<TurnInfo> compute_coiling_turns_for_chain_with_params(
    const RDKit::ROMol& polymer, int increment, int turn_size, int chain_size,
    int first_chain_size)
{

    std::vector<TurnInfo> turns;

    unsigned int residues_accounted = 0;
    auto total_monomers = polymer.getNumAtoms();
    unsigned int current_chain_size = first_chain_size;
    bool turn_downward = true;
    while (residues_accounted < total_monomers) {
        if (turn_size < 0) {
            // invalid parameters, return empty vector to indicate failure
            return {};
        }
        unsigned int turn_position = residues_accounted + current_chain_size;
        if (turn_position >= total_monomers) {
            // no more turns needed, we're done
            break;
        }
        turns.push_back({turn_position, static_cast<unsigned int>(turn_size),
                         turn_downward});
        residues_accounted += current_chain_size + turn_size;
        current_chain_size = chain_size;
        turn_size += increment;
        turn_downward = !turn_downward;
    }
    return turns;
}

float get_position_in_coils_float(unsigned int monomer_idx,
                                  const RDKit::ROMol& polymer,
                                  const std::vector<TurnInfo>& turns)
{
    int section_number = 0;
    int begin_of_section = 0;
    int section_size = 0;
    int last_segment_size = 0;

    for (unsigned int turn_idx = 0; turn_idx < turns.size(); ++turn_idx) {
        const auto& turn = turns[turn_idx];
        section_number = turn_idx * 2;
        last_segment_size = turn.position - begin_of_section;
        if (monomer_idx < turn.position + turn.size) {
            if (monomer_idx >= turn.position) {
                // we're in the turn
                section_number += 1;
                begin_of_section = turn.position;
                section_size = turn.size;
            } else {
                // we're in the segment before the turn
                section_size = last_segment_size;
            }

            return section_number +
                   static_cast<double>(monomer_idx + 1 - begin_of_section) /
                       (section_size + 1);
            break;
        }
        begin_of_section = turn.position + turn.size;
    }
    // we are in the last segment after the last turn. This segment might not be
    // complete, so the best way to approximate the position score is to assume
    // that the segment has the same size as the last complete segment
    section_size = last_segment_size;

    section_number = turns.size() * 2;

    return section_number +
           static_cast<double>(monomer_idx - begin_of_section) /
               (section_size + 1);
}

float score_coiling_layout(const RDKit::ROMol& polymer,
                           const std::vector<TurnInfo>& turns,
                           const std::vector<CustomBondInfo>& custom_bonds)
{
    if (custom_bonds.empty()) {
        return 0.0;
    }

    float DEVIATION_PENALTY =
        0.1; // penalty factor for deviation from ideal turn size

    float return_score = 0.0;
    for (auto turn : turns) {
        return_score += DEVIATION_PENALTY * std::abs(turn.y_size_difference);
    }
    for (auto bond : custom_bonds) {
        auto pos_i =
            get_position_in_coils_float(bond.monomer_i, polymer, turns);
        auto pos_j =
            get_position_in_coils_float(bond.monomer_j, polymer, turns);

        // ideal distance between bonded monomers in the coiling layout is 4
        // (one full coil), so score is based on how close the actual distance
        // is to 4. We use a quadratic function to penalize larger deviations
        // more
        double distance = std::abs(pos_j - pos_i);
        double score = (distance - 4.0) * (distance - 4.0);
        return_score += score;
    }
    return return_score / custom_bonds.size();
}

/**
 * Computes turn information for coiling chain layout.
 *
 * Analyzes CUSTOM_BONDs to determine if the polymer can be laid out as a
 * coiling pattern where chains alternate between top and bottom positions
 * (y = 0, -1, +1, -2, +2, ...), with each turn progressively longer to
 * approximate an ellipse.
 *
 * @param polymer The polymer molecule to analyze
 * @return Vector of turn information if coiling layout is feasible, empty
 *         vector otherwise
 */
static std::vector<TurnInfo>
compute_coiling_turns_for_chain(const RDKit::ROMol& polymer)
{
    // Extract all CUSTOM_BONDs
    auto custom_bonds = extract_custom_bonds(polymer);

    // Need CUSTOM_BONDs to determine pattern
    if (custom_bonds.empty()) {
        return {};
    }

    // custom bonds are ordered by monomer_i. If monomer_j are not listed
    // in ascending order, the bonds cross and the coiling layout is not
    // possible
    for (size_t idx = 1; idx < custom_bonds.size(); ++idx) {
        if (custom_bonds[idx].monomer_j < custom_bonds[idx - 1].monomer_j) {
            return {};
        }
    }

    // now we need to figure out a layout that puts connected monomers in
    // "corresponding" positions in the layers (i.e. bound monomers need to be
    // separated by exactly a full coil, or as close as possible to that
    // position). This way the bond will not cross with the coiling backbone

    // score different layouts to see if any of them would put the bonds in good
    // position. Variables are: chain_size, first_chain_size (which can be
    // smaller than the others, to allow for more flexibility in the layout),
    // first_turn_size and increment (+1 or -1). Each chain after the first one
    // will have the same size, and each turn will be increment bigger (or
    // smaller if increment is negative) than the previous one. Depending on the
    // sign of increment the spiral will expand from the center or become
    // progressively smaller (the equivalent of reversing the sequence). The
    // best layout will be the one that minimizes the distance between connected
    // monomers and their ideal positions in the layers.
    std::vector<std::vector<TurnInfo>> possible_layouts;
    for (auto increment : {2, -2}) {
        for (auto turn_size = 0; turn_size < 3; ++turn_size) {
            for (auto chain_size = 1; chain_size < 8; ++chain_size) {
                for (auto first_chain_size = 1; first_chain_size <= chain_size;
                     ++first_chain_size) {
                    auto first_turn_size = turn_size;
                    // if increment is negative, turn_size is actually the size
                    // of the last turn. Compute a plausible first turn size
                    // based on the number of residues and the chain sizes
                    if (increment < 0) {
                        int placed_residues = first_chain_size;
                        int current_turn_size = first_turn_size;
                        while (static_cast<unsigned int>(placed_residues) <
                               polymer.getNumAtoms()) {
                            placed_residues += current_turn_size + chain_size;
                            current_turn_size += abs(increment);
                        }
                        first_turn_size = current_turn_size - abs(increment);
                    }
                    auto turns = compute_coiling_turns_for_chain_with_params(
                        polymer, increment, first_turn_size, chain_size,
                        first_chain_size);
                    if (!turns.empty()) {
                        possible_layouts.push_back(turns);
                    }
                }
            }
        }
    }
    // return the best_scoring layout, or empty vector if no valid layouts found
    if (possible_layouts.empty()) {
        return {};
    }

    std::vector<std::pair<std::vector<TurnInfo>, float>> layouts_with_scores;
    std::transform(possible_layouts.begin(), possible_layouts.end(),
                   std::back_inserter(layouts_with_scores),
                   [&polymer, &custom_bonds](auto layout) {
                       return std::make_pair(
                           layout,
                           score_coiling_layout(polymer, layout, custom_bonds));
                   });
    auto min_layout_and_score =
        std::min_element(layouts_with_scores.begin(), layouts_with_scores.end(),
                         [](auto layout1, auto layout2) {
                             return layout1.second < layout2.second;
                         });
    return min_layout_and_score->first;
}

/**
 * Lays out a chain with explicit turn positions. Used for both snaking and
 * coiling layouts
 *
 * @param polymer The polymer molecule to lay out
 * @param placed_monomers_idcs Set to track which monomers have been placed
 * @param turn_positions Vector of turns specifying where turns occur, their
 * sizes and the direction to turn towards. Each turn position marks the END of
 * a segment. Turn size specifies how many residues are consumed in the turn.
 * For example, [{10, 0}, {20, 2}] means:
 *   - monomers 0-9 in first row (instant turn at 10)
 *   - monomers 10-19 in second row (turn consumes residues 20-21)
 *   - monomers 22-end in third row
 */
static void
lay_out_chain_with_turns(RDKit::ROMol& polymer,
                         std::unordered_set<int>& placed_monomers_idcs,
                         const std::vector<TurnInfo>& turn_positions)
{
    auto& conformer = polymer.getConformer();
    RDGeom::Point3D chain_start_pos(0, 0, 0);
    ChainDirection chain_dir = ChainDirection::LTR;
    auto total_monomers = polymer.getNumAtoms();

    unsigned int segment_start = 0;

    // Lay out each segment between turns
    for (size_t turn_idx = 0; turn_idx <= turn_positions.size(); ++turn_idx) {
        // Determine segment end: either the next turn position or the end of
        // chain
        unsigned int segment_end = (turn_idx < turn_positions.size())
                                       ? turn_positions[turn_idx].position
                                       : total_monomers;

        // Get turn size for this turn (0 if we're at the end)
        unsigned int turn_size = (turn_idx < turn_positions.size())
                                     ? turn_positions[turn_idx].size
                                     : 0;

        if (segment_start >= segment_end) {
            continue; // Skip empty segments
        }

        // Lay out this segment
        lay_out_chain(polymer, polymer.getAtomWithIdx(segment_start),
                      placed_monomers_idcs, chain_start_pos, chain_dir,
                      BranchDirection::UP,
                      [segment_end](const RDKit::Atom* monomer) {
                          return monomer->getIdx() == segment_end;
                      });

        // Lay out turn monomers (if any) and prepare for next segment

        RDGeom::Point3D turn_start = conformer.getAtomPos(segment_end - 1);

        // Lay out turn monomers following polygon arc
        chain_start_pos =
            lay_out_turn(polymer, conformer, placed_monomers_idcs, turn_start,
                         segment_end, turn_size, total_monomers, chain_dir,
                         turn_positions[turn_idx].downward,
                         turn_positions[turn_idx].y_size_difference);

        // Reverse direction for next segment
        chain_dir = (chain_dir == ChainDirection::LTR) ? ChainDirection::RTL
                                                       : ChainDirection::LTR;

        // Skip turn_size residues before starting the next segment
        segment_start = segment_end + turn_size;
    }
}

// TODO use Kevin's logic from
// https://github.com/schrodinger/sketcher/pull/233/changes#diff-0c7dcb7c9a99f114e01b6b000ed280af1f5ad48bbd6b0c0b83b2eda515214b55R60-R71
static bool is_backbone_bond(const RDKit::Bond* bond)
{
    if (!bond->hasProp(CUSTOM_BOND)) {
        // No custom bond property, this is a regular backbone bond
        return true;
    }

    std::string custom_bond;
    bond->getPropIfPresent(CUSTOM_BOND, custom_bond);
    return (custom_bond == "R2-R1");
}

/**
 * Helper function to check if a polymer can be considered as a single chain. If
 * there are any rings, they must be closed by non-backbone connections (e.g. a
 * turn closed by a disulfide bond), so that the main backbone can still be
 * laid out as a single chain with turns.
 */
static bool is_a_chain(const RDKit::ROMol& polymer)
{
    // Check that the polymer is a single chain (if there are rings they can be
    // closed only by non-backbone connections). We check this by counting
    // backbone bonds to monomers: exactly two monomers should have only one
    // backbone bond (the two ends), and the rest should have exactly two.

    int end_monomers = 0;
    for (const auto& atom : polymer.atoms()) {
        int num_backbone_bonds = 0;

        for (const auto& neighbor : polymer.atomNeighbors(atom)) {
            if (is_backbone_bond(polymer.getBondBetweenAtoms(
                    atom->getIdx(), neighbor->getIdx()))) {
                num_backbone_bonds++;
            }
        }
        if (num_backbone_bonds > 2) {
            return false;
        }
        if (num_backbone_bonds == 1) {
            ++end_monomers;
        }
    }
    return (end_monomers == 2);
}

/**
 * Attempts to lay out a polymer with cycles as a snaking chain. This is only
 * possible if the cycles are closed by non-backbone connections (e.g. a turn
 * closed by a disulfide bond), and if the turn constraints from multiple
 * CUSTOM_BONDS are compatible. If successful, the polymer will be laid out in
 * a snaking pattern, placed_monomers_idcs will be updated for every atom, and
 * the function will return true. If not successful, the function will return
 * false, coordinates will not be changed and placed_monomers_idcs will not be
 * updated
 */
static bool maybe_lay_out_cyclic_polymer_as_snaking_chain(
    RDKit::ROMol& polymer, std::unordered_set<int>& placed_monomers_idcs)
{
    auto turns = compute_turn_positions_for_chain(polymer);

    if (turns.size() == 0) {
        return false;
    }
    std::vector<TurnInfo> turn_positions;
    for (const auto& turn : turns) {
        // Calculate distance from the constraint range
        unsigned int distance = turn.max_pos - turn.min_pos;

        // Determine turn size based on distance, so that the CUSTOM_BOND is
        // laid out vertically
        unsigned int turn_size;
        if (distance < 4) {
            turn_size = 0; // Very small distance: instant turn
        } else if (distance % 2 == 1) {
            turn_size = 2; // Odd distance
        } else {
            turn_size = 1; // Even distance
        }

        unsigned int turn_pos =
            (turn.min_pos + 1 + (turn.max_pos - turn_size - turn.min_pos) / 2);

        turn_positions.push_back({turn_pos, turn_size});
    }
    lay_out_chain_with_turns(polymer, placed_monomers_idcs, turn_positions);
    return true;
}

/**
 * modifies the turn information by changing the size of the turn at turn_i by
 * size_change, and adjusting the positions of subsequent turns accordingly.
 * This is used to optimize the coiling layout by trying different turn sizes
 * and seeing how they affect the overall layout score.
 */
std::vector<TurnInfo> change_turn_size(const std::vector<TurnInfo>& turns,
                                       int turn_i, int size_change)
{
    std::vector<TurnInfo> new_turns = turns;
    new_turns[turn_i].size += size_change;
    new_turns[turn_i].y_size_difference -= size_change;

    for (unsigned int i = turn_i + 1; i < new_turns.size(); ++i) {
        new_turns[i].position += size_change;
    }
    return new_turns;
}

/**
 * Attempts to optimize the coiling layout by adjusting turn sizes to better
 * align the bonds with the ideal positions in the coil. The function scores
 * the initial layout and then tries increasing or decreasing the size of each
 * turn to see if it improves the score. If a better layout is found, the turns
 * are updated accordingly. The function returns true if the final layout is
 * acceptable based on a score threshold, or false if no acceptable layout could
 * be found (e.g. due to incompatible constraints from CUSTOM_BONDs).
 */
static bool optimize_coiling_layout(RDKit::ROMol& polymer,
                                    std::vector<TurnInfo>& turns)
{
    // Score the layout to evaluate bond alignment
    auto custom_bonds = extract_custom_bonds(polymer);
    auto score = score_coiling_layout(polymer, turns, custom_bonds);

    const float GOOD_SCORE = 0.1; // Threshold for acceptable layout
    if (score < GOOD_SCORE) {
        return true; // Perfect layout, no optimization needed
    }
    auto best_layout = turns;
    float best_score = score;
    // Try optimizing by adjusting turn sizes (e.g. increasing or
    // decreasing turn sizes to better alignment)

    for (unsigned int turn_i = 0; turn_i < best_layout.size(); ++turn_i) {
        for (auto change : {-1, 1}) {
            if ((turn_i > 0 && best_layout[turn_i].size + change < 2) ||
                (best_layout[turn_i].size == 0 && change < 0)) {
                continue; // use very small turns only at the beginning of the
                          // coil, to avoid clashes and avoid negative turn
                          // sizes
            }
            auto new_layout = change_turn_size(best_layout, turn_i, change);
            auto score =
                score_coiling_layout(polymer, new_layout, custom_bonds);
            if (score < best_score) {
                best_score = score;
                best_layout = new_layout;
                break;
            }
        }
    }
    turns = best_layout;
    // Threshold for acceptable layout after optimization
    const float ACCEPTABLE_SCORE = 1.5;
    return best_score < ACCEPTABLE_SCORE;
}

/**
 * Attempts to lay out a polymer with cycles as a coiling chain. This is only
 * possible if the cycles are closed by non-backbone connections (e.g. a turn
 * closed by a disulfide bond), and if the pattern of CUSTOM_BOND constraints is
 * compatible with a coiling layout. If successful, the polymer will be laid out
 * in a coiling pattern, placed_monomers_idcs will be updated for every atom,
 * and the function will return true. If not successful, the function will
 * return false, coordinates will not be changed and placed_monomers_idcs will
 * not be updated
 */
static bool maybe_lay_out_cyclic_polymer_as_coiling_chain(
    RDKit::ROMol& polymer, std::unordered_set<int>& placed_monomers_idcs)
{
    // Analyze CUSTOM_BONDs to determine if pattern fits coiling layout
    auto turns = compute_coiling_turns_for_chain(polymer);

    if (turns.empty()) {
        return false; // Pattern doesn't fit coiling layout
    }
    if (!optimize_coiling_layout(polymer, turns)) {
        return false; // Optimization failed, likely due to incompatible
                      // constraints
    }
    // Perform the layout
    lay_out_chain_with_turns(polymer, placed_monomers_idcs, turns);
    return true;
}

static void
lay_out_cyclic_polymer(RDKit::ROMol& polymer,
                       std::unordered_set<int>& placed_monomers_idcs)
{

    if (is_a_chain(polymer)) {
        // try first to lay out with a snaking layout
        if (maybe_lay_out_cyclic_polymer_as_snaking_chain(
                polymer, placed_monomers_idcs)) {
            // laid out successfully, nothing else to do
            return;
        }
        // if the first attempt fails, try laying out as a coiling chain
        if (maybe_lay_out_cyclic_polymer_as_coiling_chain(
                polymer, placed_monomers_idcs)) {
            // laid out successfully, nothing else to do
            return;
        }
    }

    // If the approach above fails, fall back to laying out the rings first
    // and then extending chains from there
    generate_coordinates_for_cycles(polymer, placed_monomers_idcs);

    for (auto ring_atoms : polymer.getRingInfo()->atomRings()) {
        for (auto monomer_idx : ring_atoms) {
            auto monomer = polymer.getAtomWithIdx(monomer_idx);
            // starting from each atom bound to a ring, lay out chains
            // horizontally to the left and right, pointing away from the ring
            for (auto neighbor : polymer.atomNeighbors(monomer)) {
                if (polymer.getRingInfo()->numAtomRings(neighbor->getIdx()) >
                    0) {
                    continue;
                }

                auto neighbor_pos =
                    polymer.getConformer().getAtomPos(neighbor->getIdx());

                auto dir = neighbor_pos -
                           polymer.getConformer().getAtomPos(monomer_idx);

                ChainDirection chain_dir =
                    dir.x < 0 ? ChainDirection::RTL : ChainDirection::LTR;
                lay_out_chain(polymer, neighbor, placed_monomers_idcs,
                              neighbor_pos, chain_dir);
            }
        }
    }
}

/**
 * Finds the "hairpin turn" in `polymer` i.e. a monomer ring with exactly one
 * zero order bond
 */
static std::vector<const RDKit::Atom*>
find_hairpin_turn(const RDKit::ROMol& polymer)
{
    // Determine the hairpin turn by finding a "ring" with only one zero order
    // bond.
    compute_full_ring_info(polymer);
    auto monomer_rings = polymer.getRingInfo()->atomRings();
    auto bond_rings = polymer.getRingInfo()->bondRings();
    std::vector<int> hairpin_turn_idxs{};
    for (size_t i = 0; i < bond_rings.size(); i++) {
        auto bond_ring = bond_rings[i];
        int zob_count = std::count_if(
            bond_ring.begin(), bond_ring.end(), [polymer](int bond_idx) {
                return polymer.getBondWithIdx(bond_idx)->getBondType() ==
                       RDKit::Bond::BondType::ZERO;
            });
        if (zob_count == 1) {
            hairpin_turn_idxs = monomer_rings[i];
            break;
        }
    }
    if (hairpin_turn_idxs.size() == 0) {
        throw std::runtime_error("Error generating coordinates for hairpin "
                                 "style polymer: no ring found "
                                 "with exactly one zero order bond.");
    }

    // collect non-branch monomers in the hairpin turn and orient the vector so
    // the start of the turn (top end of the hairpin turn) is at the front of
    // the vector
    std::vector<const RDKit::Atom*> turn_monomers{};
    auto turn_ends = std::make_pair<int, int>(-1, -1);
    for (unsigned int monomer_idx : hairpin_turn_idxs) {
        auto monomer = polymer.getAtomWithIdx(monomer_idx);
        if (monomer->getProp<bool>(BRANCH_MONOMER)) {
            continue;
        }
        turn_monomers.push_back(monomer);
        if (polymer.getRingInfo()->numAtomRings(monomer_idx) > 1) {
            if (turn_ends.first == -1) {
                turn_ends.first = turn_monomers.size() - 1;
            } else {
                turn_ends.second = turn_monomers.size() - 1;
                if (turn_monomers[turn_ends.first]->getIdx() > monomer_idx) {
                    std::swap(turn_ends.first, turn_ends.second);
                }
            }
        }
    }

    // turn_ends are the indices, in turn_monomers, of the monomers at the top
    // and bottom end of the hairpin turn respectively - so they're either
    // adjacent or at the ends of the vector. We want to orient the vector so
    // the monomer at turn_ends.first is at the front and the monomer at
    // turn_ends.second is at the back.
    auto turn_start = turn_ends.first;
    if (turn_ends.second == turn_ends.first + 1 ||
        turn_ends.second < turn_ends.first - 1) {
        std::reverse(turn_monomers.begin(), turn_monomers.end());
        turn_start = turn_monomers.size() - 1 - turn_start;
    }
    std::rotate(turn_monomers.begin(), turn_monomers.begin() + turn_start,
                turn_monomers.end());

    return turn_monomers;
}

/**
 * Lays out a "hairpin" style polymer i.e. a polymer with one or more hydrogen
 * bond connections with itself.
 */
static void
layout_hairpin_polymer(RDKit::ROMol& polymer,
                       std::unordered_set<int>& placed_monomers_idcs)
{
    auto hairpin_turn = find_hairpin_turn(polymer);
    auto& conformer = polymer.getConformer();

    // layout the hairpin turn
    // adding extra dummy vertices to the n-gon so the hairpin turn monomers are
    // laid out in a crescent from 2 to n-1 (those numbers are for 2 extra
    // vertices). We add more extra vertices for turns with less than 8 monomers
    // to increase their radius enough so that the chains coming off the hairpin
    // turn have enough vertical distance between them to comfortably display
    // the bonds between branch monomers in those chains.
    auto dummy_vertex_count = hairpin_turn.size() < 8 ? 4 : 2;
    auto cycle_coords = get_coords_for_ngon(
        hairpin_turn.size() + dummy_vertex_count, PolygonStartSide::LEFT);
    int turn_start = (dummy_vertex_count / 2) + 1;
    int turn_end = turn_start + hairpin_turn.size() - 1;
    // lay out all the monomers in the hairpin turn except the ends - they'll
    // get laid out when we lay out the chains
    for (size_t i = 1; i < hairpin_turn.size() - 1; i++) {
        auto monomer = hairpin_turn[i];
        auto coords = cycle_coords[i + turn_start];
        place_monomer_at(conformer, monomer, coords, placed_monomers_idcs);

        for (auto neighbor : polymer.atomNeighbors(monomer)) {
            if (neighbor->getProp<bool>(BRANCH_MONOMER)) {
                auto radius = coords.length();
                auto multiplier = (radius + MONOMER_BOND_LENGTH) / radius;
                place_monomer_at(conformer, neighbor, coords * multiplier,
                                 placed_monomers_idcs);
            }
        }
    }

    // layout the chains
    lay_out_chain(polymer, hairpin_turn.front(), placed_monomers_idcs,
                  cycle_coords[turn_start], ChainDirection::RTL,
                  BranchDirection::DOWN);
    lay_out_chain(polymer, hairpin_turn.back(), placed_monomers_idcs,
                  cycle_coords[turn_end], ChainDirection::RTL,
                  BranchDirection::UP);
}

static void
lay_out_linear_polymer(RDKit::ROMol& polymer,
                       std::unordered_set<int>& placed_monomers_idcs,
                       const bool rotate = false)
{
    auto monomers = polymer.atoms();
    auto start_monomer =
        *std::find_if(monomers.begin(), monomers.end(), [](RDKit::Atom* m) {
            return !m->getProp<bool>(BRANCH_MONOMER);
        });
    auto branch_dir = rotate ? BranchDirection::UP : BranchDirection::DOWN;
    auto chain_dir = rotate ? ChainDirection::RTL : ChainDirection::LTR;
    lay_out_chain(polymer, start_monomer, placed_monomers_idcs,
                  RDGeom::Point3D(0, 0, 0), chain_dir, branch_dir);
}

/**
 * Lays out a polymer in the appropriate layout based on its monomer bonds.
 * Supports the following layouts:
 * * Linear:
 *   .__.__.__.__.__.__.
 *      |     |        |
 *      *     *        *
 *
 * * Cyclic:
 *     __
 *    /  \__.__.__.
 *   |   |__.
 *   \__/
 *
 * * Hairpin:
 *              __
 *       .__.__/  \
 *   .__.__.__    |
 *            \__/
 *
 * @param polymer
 * @param rotate This only applies for linear polymers. Rotates the direction of
 * monomer chains and branches 180° (chains go RTL and branches go UP).
 */
static void lay_out_polymer(RDKit::ROMol& polymer,
                            std::unordered_set<int>& placed_monomers_idcs,
                            const bool rotate = false)
{
    if (!polymer.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(polymer, /*res=*/nullptr, include_dative_bonds);
    }
    if (polymer.getRingInfo()->numRings() > 0) {
        lay_out_cyclic_polymer(polymer, placed_monomers_idcs);
    } else if (polymer.getNumAtoms() == polymer.getNumBonds() + 1) {
        lay_out_linear_polymer(polymer, placed_monomers_idcs, rotate);
    } else {
        layout_hairpin_polymer(polymer, placed_monomers_idcs);
    }
}

/**
 * Lays out a simple long linear polymer (with no branches) in a snaking
 * pattern. Distributes monomers as evenly as possible across rows to avoid
 * having a very short last row.
 *
 * The algorithm uses both 1-residue and 2-residue turns to find the optimal
 * layout that:
 *   1. Keeps all straight segments <= MAX_MONOMERS_PER_SNAKE
 *   2. Distributes monomers evenly across rows
 */
static void lay_out_snaked_linear_polymer(RDKit::ROMol& polymer)
{
    auto placed_monomers_idcs = std::unordered_set<int>{};
    auto total_monomers = polymer.getNumAtoms();

    unsigned int num_rows =
        (total_monomers + MAX_MONOMERS_PER_SNAKE - 1) / MAX_MONOMERS_PER_SNAKE;
    unsigned int turn_size = 1;

    // Calculate  number of rows needed if all segments MAX_MONOMERS_PER_SNAKE
    // With num_rows segments, we have (num_rows - 1) turns between
    // them. Calculate how much space is left for straight segments after
    // accounting for turns
    int segment_space = total_monomers - (num_rows - 1) * turn_size;

    // Distribute the segment space as evenly as possible across all segments.
    // Monomers that are left will be accommodated by using increased turn sizes
    // for as many rows as needed.
    unsigned int segment_length = segment_space / num_rows;
    unsigned int num_of_longer_turns = segment_space % num_rows;
    unsigned int current_pos = 0;

    std::vector<TurnInfo> turn_positions;
    for (unsigned int row = 0; row < num_rows; ++row) {
        current_pos += segment_length;
        unsigned int turn_size = (row < num_of_longer_turns) ? 2 : 1;
        // turns always take space like a 1 residue turn, so they can be aligned
        int y_size_difference =
            1 - turn_size; // -1 for 2-residue turns, 0 for 1-residue turns

        // Add a turn after this segment (except for the last segment)
        if (row < num_rows - 1) {
            turn_positions.push_back(
                {current_pos, turn_size, true, y_size_difference});
            // Advance position past the turn residues
            current_pos += turn_size;
        }
    }
    // Use the calculated turn positions to lay out the snaking chain
    lay_out_chain_with_turns(polymer, placed_monomers_idcs, turn_positions);
}

/**
 * Sorts polymers in connection order i.e. a polymer is followed by all its
 * neighbors contiguously that are then followed by their neighbors contiguously
 * and so on.
 * @return the sorted polymers vector and a map of child polymer to parent
 * polymer
 */
std::pair<std::vector<RDKit::ROMOL_SPTR>,
          std::map<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>>
sort_polymers_by_connectivity(const std::vector<RDKit::ROMOL_SPTR>& polymers)
{
    std::map<unsigned int, RDKit::ROMOL_SPTR> polymer_for_monomer_idx;
    for (auto polymer : polymers) {
        for (auto monomer : polymer->atoms()) {
            auto monomer_idx = monomer->getProp<unsigned int>(ORIGINAL_INDEX);
            polymer_for_monomer_idx[monomer_idx] = polymer;
        }
    }

    std::vector<RDKit::ROMOL_SPTR> sorted_polymers;
    std::set<std::string> visited_polymers_ids;
    std::map<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR> parent_polymer;
    for (auto polymer : polymers) {
        std::queue<RDKit::ROMOL_SPTR> polymer_queue;
        for (polymer_queue.push(polymer); !polymer_queue.empty();
             polymer_queue.pop()) {
            auto polymer = polymer_queue.front();
            auto polymer_id = polymer->getProp<std::string>(POLYMER_ID);
            if (visited_polymers_ids.contains(polymer_id)) {
                continue;
            }
            visited_polymers_ids.insert(polymer_id);
            sorted_polymers.push_back(polymer);

            for (auto monomer : polymer->atoms()) {
                if (!monomer->hasProp(BOND_TO)) {
                    continue;
                }
                std::vector<int> neighbor_monomer_idcs;
                monomer->getPropIfPresent<std::vector<int>>(
                    BOND_TO, neighbor_monomer_idcs);
                for (auto neighbor_monomer_idx : neighbor_monomer_idcs) {
                    auto neighbor_polymer =
                        polymer_for_monomer_idx[neighbor_monomer_idx];
                    if (!visited_polymers_ids.contains(
                            neighbor_polymer->getProp<std::string>(
                                POLYMER_ID))) {
                        parent_polymer[neighbor_polymer] = polymer;
                    }
                    polymer_queue.push(neighbor_polymer);
                }
            }
        }
    }

    return std::make_pair(sorted_polymers, parent_polymer);
}

static std::pair<std::vector<RDKit::ROMOL_SPTR>,
                 std::map<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>>
break_into_polymers(const RDKit::ROMol& monomer_mol)
{
    for (auto atom : monomer_mol.atoms()) {
        atom->setProp(ORIGINAL_INDEX, atom->getIdx());
    }

    for (auto bond : monomer_mol.bonds()) {
        auto beginMonomer = bond->getBeginAtom();
        auto endMonomer = bond->getEndAtom();
        if (get_polymer_id(beginMonomer) == get_polymer_id(endMonomer)) {
            continue;
        }
        // store the connected indices as comma separated values in the
        // BOND_TO prop
        std::vector<int> begin_monomer_bond_to;
        if (beginMonomer->hasProp(BOND_TO)) {
            begin_monomer_bond_to =
                beginMonomer->getProp<std::vector<int>>(BOND_TO);
        }
        begin_monomer_bond_to.push_back(endMonomer->getIdx());
        beginMonomer->setProp(BOND_TO, begin_monomer_bond_to);
        std::vector<int> end_monomer_bond_to;
        if (endMonomer->hasProp(BOND_TO)) {
            end_monomer_bond_to =
                endMonomer->getProp<std::vector<int>>(BOND_TO);
        }
        end_monomer_bond_to.push_back(beginMonomer->getIdx());
        endMonomer->setProp(BOND_TO, end_monomer_bond_to);
    }

    std::vector<RDKit::ROMOL_SPTR> polymers;
    for (auto polymer_id : get_polymer_ids(monomer_mol)) {
        auto polymer = extract_helm_polymers(monomer_mol, {polymer_id});
        polymer->setProp(POLYMER_ID, polymer_id);
        polymer->addConformer(new RDKit::Conformer(polymer->getNumAtoms()));
        polymers.push_back(polymer);
    }

    return sort_polymers_by_connectivity(polymers);
}

typedef std::vector<std::pair<unsigned int, unsigned int>> BOND_IDX_VEC;
static BOND_IDX_VEC get_bonds_between_polymers(const RDKit::ROMol& from,
                                               const RDKit::ROMol& to)
{
    std::map<unsigned int, RDKit::Atom*> end_monomers;
    for (auto monomer : to.atoms()) {
        end_monomers[monomer->getProp<unsigned int>(ORIGINAL_INDEX)] = monomer;
    }

    BOND_IDX_VEC bonds{};
    for (auto monomer : from.atoms()) {
        std::vector<int> bonded_monomers_indices;
        monomer->getPropIfPresent<std::vector<int>>(BOND_TO,
                                                    bonded_monomers_indices);
        if (bonded_monomers_indices.empty()) {
            continue;
        }
        for (auto bonded_idx : bonded_monomers_indices) {
            if (end_monomers.contains(bonded_idx)) {
                bonds.push_back(
                    {monomer->getIdx(), end_monomers[bonded_idx]->getIdx()});
            }
        }
    }
    return bonds;
}
// Assign unit IDs to ring systems. Rings that share at least one atom get the
// same unit ID. This ensures that fused ring systems are treated as a single
// topological unit.
static void assign_ring_unit_ids(const RDKit::ROMol& mol,
                                 std::unordered_map<int, std::string>& unit_ids)
{
    if (!mol.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(mol, /*res=*/nullptr, include_dative_bonds);
    }

    // Assign each ring its own unique ID
    // Each ring initially gets an ID based on its first atom index
    for (auto cycle : mol.getRingInfo()->atomRings()) {
        std::string ring_id = "ring_" + std::to_string(cycle[0]);
        for (auto atom_idx : cycle) {
            unit_ids[atom_idx] = ring_id;
        }
    }

    // Step 2: Iteratively merge rings that share atoms
    // We need to make sure that if A shares with B and B shares with C, all
    // three get the same ID
    auto rings = mol.getRingInfo()->atomRings();
    bool changed = true;

    while (changed) {
        changed = false;

        // Check each pair of rings
        for (size_t i = 0; i < rings.size(); ++i) {
            for (size_t j = i + 1; j < rings.size(); ++j) {
                // Get current IDs for these rings (check the first atom of
                // each)
                std::string id_i = unit_ids[rings[i][0]];
                std::string id_j = unit_ids[rings[j][0]];

                // Skip if already merged
                if (id_i == id_j) {
                    continue;
                }

                // Check if rings share any atoms
                bool share_atom = false;
                for (auto atom_i : rings[i]) {
                    for (auto atom_j : rings[j]) {
                        if (atom_i == atom_j) {
                            share_atom = true;
                            break;
                        }
                    }
                    if (share_atom)
                        break;
                }

                if (share_atom) {
                    // Merge: pick lexicographically smaller ID as canonical
                    // This ensures a stable, deterministic choice
                    std::string canonical = std::min(id_i, id_j);
                    std::string to_replace = std::max(id_i, id_j);

                    // Copy canonical tag to all atoms with to_replace tag
                    for (auto& [atom_idx, ring_id] : unit_ids) {
                        if (ring_id == to_replace) {
                            ring_id = canonical;
                            changed = true;
                        }
                    }
                }
            }
        }
    }
}

static std::pair<std::vector<std::pair<int, int>>, std::unordered_map<int, int>>
find_chains(const RDKit::ROMol& mol, int start_atom_idx, int ring_atom_id,
            std::unordered_map<int, std::string>& unit_ids)
{
    // BFS to find chains off a ring atom
    std::unordered_map<int, bool> visited;
    std::unordered_map<int, int> distance;
    std::unordered_map<int, int> predecessor;
    std::queue<int> atom_queue;

    atom_queue.push(start_atom_idx);
    visited[start_atom_idx] = true;
    distance[start_atom_idx] = 0;
    predecessor[start_atom_idx] = -1;

    while (!atom_queue.empty()) {
        int atom_idx = atom_queue.front();
        atom_queue.pop();
        auto atom = mol.getAtomWithIdx(atom_idx);
        for (auto nbr : mol.atomNeighbors(atom)) {
            int nbr_idx = nbr->getIdx();
            if (unit_ids[nbr_idx].empty() && !visited[nbr_idx]) {
                atom_queue.push(nbr_idx);
                predecessor[nbr_idx] = atom_idx;
                distance[nbr_idx] = distance[atom_idx] + 1;
                visited[nbr_idx] = true;
            }
        }
    }
    std::vector<std::pair<int, int>> sorted_distances(distance.begin(),
                                                      distance.end());
    std::sort(sorted_distances.begin(), sorted_distances.end(),
              [](auto& a, auto& b) { return a.second > b.second; });
    return {sorted_distances, predecessor};
}

// Assign unit IDs to chains connected to rings. Chains directly connected to
// rings get the same unit ID as the ring, branches off those chains get new
// unit IDs. Takes a map of atom idx to unit ID string to fill it.
static void
assign_chain_unit_ids(const RDKit::ROMol& mol,
                      std::unordered_map<int, std::string>& unit_ids)
{
    for (auto cycle : mol.getRingInfo()->atomRings()) {
        for (auto ring_atom_id : cycle) {
            for (auto neighbor :
                 mol.atomNeighbors(mol.getAtomWithIdx(ring_atom_id))) {
                if (unit_ids[neighbor->getIdx()].empty()) {
                    auto start_atom_idx = neighbor->getIdx();
                    auto [sorted_distances, predecessor] = find_chains(
                        mol, start_atom_idx, ring_atom_id, unit_ids);

                    // Assign chain IDs, longest gets same as ring, branches
                    // get new IDs

                    bool found_atom = true;
                    int chain_n = 1;
                    std::string chain_id = unit_ids[ring_atom_id];

                    while (found_atom) {
                        found_atom = false;
                        for (auto& dp : sorted_distances) {
                            unsigned int farthest_atom = dp.first;
                            while (farthest_atom != start_atom_idx) {
                                if (!unit_ids[farthest_atom].empty()) {
                                    break;
                                }
                                unit_ids[farthest_atom] = chain_id;
                                farthest_atom = predecessor[farthest_atom];
                                found_atom = true;
                            }
                            if (unit_ids[farthest_atom].empty()) {
                                unit_ids[farthest_atom] = chain_id;
                            }
                            chain_id = "chain_" + std::to_string(chain_n++);
                        }
                    }
                }
            }
        }
    }
}

// Set the chain IDs on atoms. Takes a map of atom idx to unit ID string that
// contains the string to set for each atom.
static void
set_atom_chain_ids(const RDKit::ROMol& mol,
                   const std::unordered_map<int, std::string>& unit_ids)
{
    if (unit_ids.size() != mol.getNumAtoms()) {
        std::cerr << "set_atom_chain_ids:: Warning: unit_ids size "
                  << unit_ids.size() << " does not match number of atoms "
                  << mol.getNumAtoms() << "\n";
        return;
    }
    for (auto atom : mol.atoms()) {
        if (auto res_info = static_cast<RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo())) {
            res_info->setChainId(unit_ids.at(atom->getIdx()));
        }
    }
}

/** assigns unit IDs to each atoms. These are strings that describes the
 * topology the atom is part of. Each ring will have a unique unit ID, chains
 * connected to rings will have the same unit ID as the ring, branches off those
 * chains will have new unit IDs. These IDs will then be used as if they were
 * chain IDs so that each topological unit is treated as a separate polymer
 * chain.
 * @return a pair containing a vector of the topological units and a map of
 * child unit to parent unit
 */
static std::pair<std::vector<RDKit::ROMOL_SPTR>,
                 std::map<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>>
break_into_topological_units(const RDKit::ROMol& monomer_mol)
{
    std::unordered_map<int, std::string> unit_ids(monomer_mol.getNumAtoms());
    assign_ring_unit_ids(monomer_mol, unit_ids);
    assign_chain_unit_ids(monomer_mol, unit_ids);
    set_atom_chain_ids(monomer_mol, unit_ids);
    return break_into_polymers(monomer_mol);
}

/**
 * Provides the y-coordinate offset to apply to each monomer in
 * `polymer_to_translate` to place it one `MONOMER_BOND_LENGTH` away below
 * `reference_polymer` with no overlap.
 */
static double get_y_offset_for_polymer(const RDKit::ROMol& polymer_to_translate,
                                       const RDKit::ROMol& reference_polymer)
{
    auto reference_coords = reference_polymer.getConformer().getPositions();
    auto lowest_point = std::min_element(
        reference_coords.begin(), reference_coords.end(),
        [](RDGeom::Point3D a, RDGeom::Point3D b) { return a.y < b.y; });

    auto to_translate_coords =
        polymer_to_translate.getConformer().getPositions();
    auto highest_point = std::max_element(
        to_translate_coords.begin(), to_translate_coords.end(),
        [](RDGeom::Point3D a, RDGeom::Point3D b) { return a.y < b.y; });

    return (*lowest_point).y - (*highest_point).y - MONOMER_BOND_LENGTH;
}

/**
 * Translate `polymer` so that it is positioned below `reference`  without
 * overlap, while preserving its internal geometry.
 *
 * If there are bonds between the two polymers, the first bonded monomer
 * of `polymer` is aligned vertically below the corresponding monomer
 * in `reference`, so that the inter-polymer bond is vertical.
 *
 * If there are no inter-polymer bonds, the polymer is simply shifted
 * downward by a fixed spacing.
 *
 * This function performs translation only (no rotation or reflection).
 * It may be called multiple times after flips to restore correct placement.
 * @param polymer the polymer to be placed
 * @param reference the reference polymer used for placement
 * @param polymer_bonds list of inter-polymer bonds. The start atom of each bond
 * is in `polymer`, the end atom is in `reference`
 */
static void place_polymer_below_reference(RDKit::ROMol& polymer,
                                          const RDKit::ROMol& reference,
                                          const BOND_IDX_VEC& polymer_bonds)
{
    // Translation offset to be applied to all atom positions
    auto pos_offset = RDGeom::Point3D(0, 0, 0);

    // Vertical separation needed to avoid overlap with the reference  polymer
    auto y_offset = get_y_offset_for_polymer(polymer, reference);

    if (polymer_bonds.empty()) {
        // No explicit connection: stack polymer below the reference
        // using a fixed vertical separation
        pos_offset.y = y_offset - DIST_BETWEEN_MULTIPLE_POLYMERS;
    } else {
        // Align the first bonded monomer so that the inter-polymer bond
        // is vertical and points downward
        auto& monomer_coords =
            polymer.getConformer().getAtomPos(polymer_bonds[0].first);
        auto& placed_monomer_coords =
            reference.getConformer().getAtomPos(polymer_bonds[0].second);

        // Target position for the bonded monomer after translation
        auto translated_monomer_coords = RDGeom::Point3D(
            placed_monomer_coords.x, monomer_coords.y + y_offset,
            placed_monomer_coords.z);

        pos_offset = translated_monomer_coords - monomer_coords;
    }

    // Apply the translation uniformly to all atoms
    for (auto& pos : polymer.getConformer().getPositions()) {
        pos += pos_offset;
    }
}

/**
 * Determine whether `polymer` needs to be flipped vertically (along the Y-axis)
 * to avoid geometric clashes between inter-polymer bonds (polymer<-> reference)
 * Special handling is applied for adjacent bonds sharing an atom, since they
 * represent a chemically valid continuation and require a more local proximity
 * check.
 *
 * Assumes that the polymer has already been placed relative to the reference.
 *
 * Returns true if a vertical flip is required to resolve overlaps.
 */
static bool
needs_vertical_flip_due_to_clashes(const RDKit::ROMol& polymer,
                                   const RDKit::ROMol& reference,
                                   const BOND_IDX_VEC& polymer_bonds)
{
    for (auto bond_between_polymers : polymer_bonds) {
        // Geometry of the inter-polymer bond
        auto interpolymer_bond_start_pos =
            polymer.getConformer().getAtomPos(bond_between_polymers.first);
        auto interpolymer_bond_end_pos =
            reference.getConformer().getAtomPos(bond_between_polymers.second);

        // Check against all bonds of the polymer
        for (auto bond_in_polymer : polymer.bonds()) {
            auto intrapolymer_bond_start_pos =
                polymer.getConformer().getAtomPos(
                    bond_in_polymer->getBeginAtomIdx());
            auto intrapolymer_bond_end_pos = polymer.getConformer().getAtomPos(
                bond_in_polymer->getEndAtomIdx());

            // Case 1: bonds sharing an atom – use adjacent bond check
            if (bond_between_polymers.first ==
                    bond_in_polymer->getBeginAtomIdx() ||
                bond_between_polymers.first ==
                    bond_in_polymer->getEndAtomIdx()) {

                auto pos1 = intrapolymer_bond_start_pos;
                auto pos2 = intrapolymer_bond_end_pos;
                if (bond_between_polymers.first ==
                    bond_in_polymer->getBeginAtomIdx()) {
                    pos1 = intrapolymer_bond_end_pos;
                    pos2 = intrapolymer_bond_start_pos;
                }

                if (adjacent_bonds_are_too_close(pos1, pos2,
                                                 interpolymer_bond_end_pos)) {
                    return true;
                }
            }
            // Case 2: non-adjacent bonds – generic proximity check
            else if (bonds_are_too_close(interpolymer_bond_start_pos,
                                         interpolymer_bond_end_pos,
                                         intrapolymer_bond_start_pos,
                                         intrapolymer_bond_end_pos)) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Determine whether `polymer` needs to be flipped horizontally (along the
 * X-axis) to resolve topological crossings between multiple inter-polymer
 * bonds.
 *
 * This uses an index-based inversion test rather than geometry and assumes:
 *
 *   - Atom indices are assumed to follow a left-to-right (LTR) layout
 *   - Atom are assumed to be listed in growing index order
 *   - Each inter-polymer bond maps a polymer atom index to a reference
 *     atom index
 *   - If the sequence of reference indices is not monotonic when sorted
 *     by polymer indices, at least one pair of connections must cross
 *
 *
 * Returns true if a horizontal flip is required to remove crossings. The flip
 * is only performed if it solves the inversion, otherwise keep the original
 * orientation even if there are fewer crossings because we assume that LTR is a
 * more readable layout.
 */
static bool needs_horizontal_flip_due_to_inversion(BOND_IDX_VEC polymer_bonds)
{
    if (polymer_bonds.size() < 2)
        return false;

    // Sort by polymer atom index (left-to-right backbone order)
    std::sort(polymer_bonds.begin(), polymer_bonds.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Flip only if the reference indices are strictly decreasing
    for (size_t i = 1; i < polymer_bonds.size(); ++i) {
        if (polymer_bonds[i].second > polymer_bonds[i - 1].second) {
            return false; // not fully reversed
        }
    }

    return true;
}

/**
 * Flip all atoms of `polymer` vertically (Y-axis) around Y = 0.
 */
static void flip_vertically(RDKit::ROMol& polymer)
{
    auto& conf = polymer.getConformer();
    for (auto& pos : conf.getPositions()) {
        pos.y = -pos.y;
    }
}

/**
 * Flip all atoms of `polymer` horizontally (X-axis) around X = 0.
 */
static void flip_horizontally(RDKit::ROMol& polymer)
{
    auto& conf = polymer.getConformer();
    for (auto& pos : conf.getPositions()) {
        pos.x = -pos.x;
    }
}

/**
 * Move `polymer_to_orient` so that it is positioned below `reference_polymer`
 * with no overlap, resolving vertical clashes and horizontal crossings.
 *
 * - First places the polymer below the reference
 * - Flips vertically if needed to resolve bond clashes
 * - Flips horizontally if needed to remove topological crossings
 * Placement is recomputed after each flip to restore correct alignment.
 * If the polymer was rotated during layout (for double stranded nucleic
 * acids), flipping steps are skipped
 */
static void orient_polymer(RDKit::ROMol& polymer_to_orient,
                           const RDKit::ROMol& reference_polymer,
                           bool polymer_was_rotated)
{
    auto polymer_bonds =
        get_bonds_between_polymers(polymer_to_orient, reference_polymer);

    // Initial placement
    place_polymer_below_reference(polymer_to_orient, reference_polymer,
                                  polymer_bonds);

    // if the polymer was rotated during layout, skip flipping steps to avoid
    // undoing the rotation
    if (polymer_was_rotated) {
        return;
    }

    // Resolve vertical clashes
    if (needs_vertical_flip_due_to_clashes(polymer_to_orient, reference_polymer,
                                           polymer_bonds)) {
        flip_vertically(polymer_to_orient);
        place_polymer_below_reference(polymer_to_orient, reference_polymer,
                                      polymer_bonds);
    }

    // Resolve horizontal crossings (index inversion)
    if (needs_horizontal_flip_due_to_inversion(polymer_bonds)) {
        flip_horizontally(polymer_to_orient);
        place_polymer_below_reference(polymer_to_orient, reference_polymer,
                                      polymer_bonds);
    }
}

/**
 * Adds a conformer to the monomer_mol copying all the coordinates from the
 * polymers. Expects each monomer in the polymers to have an ORIGINAL_INDEX
 * prop pointing to the corresponding index for that monomer in the monomer_mol.
 * Returns the id of the added conformer.
 */
static unsigned int copy_polymer_coords_to_monomer_mol(
    RDKit::ROMol& monomer_mol, const std::vector<RDKit::ROMOL_SPTR>& polymers)
{
    auto conformer = new RDKit::Conformer(monomer_mol.getNumAtoms());
    for (auto polymer : polymers) {
        for (auto monomer : polymer->atoms()) {
            conformer->setAtomPos(
                monomer->getProp<unsigned int>(ORIGINAL_INDEX),
                polymer->getConformer().getAtomPos(monomer->getIdx()));
        }
    }
    monomer_mol.addConformer(conformer);
    return conformer->getId();
}

static void remove_cxsmiles_labels(RDKit::ROMol& monomer_mol)
{
    for (auto monomer : monomer_mol.atoms()) {
        auto label = monomer->getProp<std::string>(ATOM_LABEL);
        if (label.size() > 4) {
            monomer->setProp<std::string>(ATOM_LABEL, SMILES_MONOMER_LABEL);
        }
    }
}

/**
 * Adjust CHEM polymers (polymers with a single monomer typically used to
 * connect a polymer back to itself or two polymers together) to position them
 * between the two polymers they are connecting.
 */
static void adjust_chem_polymer_coords(RDKit::ROMol& monomer_mol)
{
    auto& conformer = monomer_mol.getConformer();
    for (auto monomer : monomer_mol.atoms()) {
        if (!boost::starts_with(get_polymer_id(monomer), "CHEM")) {
            continue;
        }
        auto neighbors = monomer_mol.atomNeighbors(monomer);
        if (std::distance(neighbors.begin(), neighbors.end()) != 2) {
            continue;
        }

        auto first_neighbor = *neighbors.begin();
        auto second_neighbor = *std::next(neighbors.begin(), 1);
        if (get_polymer_id(first_neighbor) == get_polymer_id(second_neighbor)) {
            continue;
        }
        auto first_neighbor_pos =
            conformer.getAtomPos(first_neighbor->getIdx());
        auto second_neighbor_pos =
            conformer.getAtomPos(second_neighbor->getIdx());
        auto new_pos = (first_neighbor_pos + second_neighbor_pos) / 2;

        double x_offset = 0;
        if (first_neighbor_pos.x == second_neighbor_pos.x) {
            if (round(first_neighbor_pos.x) == 0) {
                x_offset = -MONOMER_BOND_LENGTH;
            } else if (first_neighbor->getDegree() == 2 &&
                       second_neighbor->getDegree() == 2) {
                x_offset = MONOMER_BOND_LENGTH;
            }
        }
        new_pos.x += x_offset;
        conformer.setAtomPos(monomer->getIdx(), new_pos);
    }
}

/**
 * Clears props that were set on monomers just to help with layout logic
 */
static void clear_layout_props(RDKit::ROMol& monomer_mol)
{
    std::vector<std::string> layout_props{BOND_TO, ORIGINAL_INDEX};
    for (auto monomer : monomer_mol.atoms()) {
        for (auto prop : layout_props) {
            monomer->clearProp(prop);
        }
    }
}

/**
 * Determines if a monomer mol contains a single linear polymer with no branches
 * or rings
 */
static bool is_single_linear_polymer(const RDKit::ROMol& monomer_mol)
{
    auto polymer_ids = get_polymer_ids(monomer_mol);
    if (polymer_ids.size() != 1) {
        return false;
    }

    compute_full_ring_info(monomer_mol);
    if (monomer_mol.getRingInfo()->numRings() > 0) {
        // reset ring info so it can get recomputed by the appropriate layout
        // method
        monomer_mol.getRingInfo()->reset();
        return false;
    }

    for (auto monomer : monomer_mol.atoms()) {
        if (monomer->getProp<bool>(BRANCH_MONOMER)) {
            return false;
        }
    }
    return true;
}

static bool are_double_stranded_nucleic_acid(RDKit::ROMOL_SPTR polymer1,
                                             RDKit::ROMOL_SPTR polymer2)
{
    // return true if the two polymers are nucleic acids and there are at least
    // two bonds between them.

    for (auto polymer : {polymer1, polymer2}) {
        if (!is_nucleic_acid(*polymer)) {
            return false;
        }
    }
    const auto polymer_bonds = get_bonds_between_polymers(*polymer1, *polymer2);
    return polymer_bonds.size() >= 2;
}

void lay_out_polymers(
    const std::vector<RDKit::ROMOL_SPTR>& polymers,
    const std::map<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>& parent_polymer)
{
    std::unordered_set<int> placed_monomers_idcs{};
    RDKit::ROMOL_SPTR last_placed_polymer = nullptr;

    // lay out the polymers in connection order so connected polymers are laid
    // out next to each other.
    for (auto polymer : polymers) {
        // If a polymer has no entry in parent_polymer, it's either the first to
        // be placed or it is not connected to any other polymer in the HELM
        // graph. In the second case, use the last placed polymer as the
        // parent/reference. This causes successive unconnected polymers to be
        // positioned relative to the previous one (typically stacked
        // vertically), instead of all being placed relative to a null parent at
        // the origin, which would make them overlap.
        RDKit::ROMOL_SPTR parent = last_placed_polymer;
        if (parent_polymer.contains(polymer)) {
            parent = parent_polymer.at(polymer);
        }
        // For double stranded nucleic acids we want to lay out the first
        // polymer normally and then rotate the other polymer 180° so the
        // strands run anti-parallel to each other
        bool rotate_polymer = (parent != nullptr) &&
                              are_double_stranded_nucleic_acid(polymer, parent);
        lay_out_polymer(*polymer, placed_monomers_idcs, rotate_polymer);
        if (parent != nullptr) {
            orient_polymer(*polymer, *parent, rotate_polymer);
        }
        last_placed_polymer = polymer;
    }
}

bool has_no_clashes(const RDKit::ROMol& monomer_mol)
{
    const auto& conformer = monomer_mol.getConformer();
    auto positions = conformer.getPositions();
    for (size_t i = 0; i < positions.size(); i++) {
        for (size_t j = i + 1; j < positions.size(); j++) {
            if (positions[i].x - positions[j].x > MONOMER_CLASH_DISTANCE ||
                positions[i].x - positions[j].x < -MONOMER_CLASH_DISTANCE ||
                positions[i].y - positions[j].y > MONOMER_CLASH_DISTANCE ||
                positions[i].y - positions[j].y < -MONOMER_CLASH_DISTANCE) {
                continue;
            }
            return false;
        }
    }
    return true;
}

bool segments_intersect(const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
                        const RDGeom::Point3D& q1, const RDGeom::Point3D& q2)
{
    segment_t seg_p = {{p1.x, p1.y}, {p2.x, p2.y}};
    segment_t seg_q = {{q1.x, q1.y}, {q2.x, q2.y}};
    std::vector<point_t> intersection_pts;
    bg::intersection(seg_p, seg_q, intersection_pts);
    return !intersection_pts.empty();
}

double compute_angle_between_adjacent_segments(const RDGeom::Point3D& a,
                                               const RDGeom::Point3D& b,
                                               const RDGeom::Point3D& c)
{
    // Vector BA = A - B
    const double bax = a.x - b.x;
    const double bay = a.y - b.y;

    // Vector BC = C - B
    const double bcx = c.x - b.x;
    const double bcy = c.y - b.y;

    // Dot and 2D cross-product magnitude
    const double dot = bax * bcx + bay * bcy;
    const double cross = bax * bcy - bay * bcx;

    return std::fabs(std::atan2(cross, dot));
}

double compute_distance_between_segments(const RDGeom::Point3D& p1,
                                         const RDGeom::Point3D& p2,
                                         const RDGeom::Point3D& q1,
                                         const RDGeom::Point3D& q2)
{
    segment_t seg_p = {{p1.x, p1.y}, {p2.x, p2.y}};
    segment_t seg_q = {{q1.x, q1.y}, {q2.x, q2.y}};
    return bg::distance(seg_p, seg_q);
}

bool adjacent_bonds_are_too_close(const RDGeom::Point3D& point1,
                                  const RDGeom::Point3D& point2,
                                  const RDGeom::Point3D& point3)
{
    return compute_angle_between_adjacent_segments(point1, point2, point3) <
           MIN_ANGLE_BETWEEN_ADJACENT_BONDS;
}

bool bonds_are_too_close(const RDGeom::Point3D& begin1_pos,
                         const RDGeom::Point3D& end1_pos,
                         const RDGeom::Point3D& begin2_pos,
                         const RDGeom::Point3D& end2_pos)
{
    auto minx1 = std::min(begin1_pos.x, end1_pos.x) - BOND_CLASH_DISTANCE;
    auto maxx1 = std::max(begin1_pos.x, end1_pos.x) + BOND_CLASH_DISTANCE;
    auto miny1 = std::min(begin1_pos.y, end1_pos.y) - BOND_CLASH_DISTANCE;
    auto maxy1 = std::max(begin1_pos.y, end1_pos.y) + BOND_CLASH_DISTANCE;
    auto minx2 = std::min(begin2_pos.x, end2_pos.x) - BOND_CLASH_DISTANCE;
    auto maxx2 = std::max(begin2_pos.x, end2_pos.x) + BOND_CLASH_DISTANCE;
    auto miny2 = std::min(begin2_pos.y, end2_pos.y) - BOND_CLASH_DISTANCE;
    auto maxy2 = std::max(begin2_pos.y, end2_pos.y) + BOND_CLASH_DISTANCE;
    // bounding box check
    if (maxx1 < minx2 || maxx2 < minx1 || maxy1 < miny2 || maxy2 < miny1) {
        return false;
    }
    return (segments_intersect(begin1_pos, end1_pos, begin2_pos, end2_pos) ||
            compute_distance_between_segments(begin1_pos, end1_pos, begin2_pos,
                                              end2_pos) < BOND_CLASH_DISTANCE);
}

bool has_no_bond_crossings(const RDKit::ROMol& monomer_mol)
{
    const auto& conformer = monomer_mol.getConformer();
    for (size_t i = 0; i < monomer_mol.getNumBonds(); i++) {
        auto bond1 = monomer_mol.getBondWithIdx(i);
        auto begin1_pos = conformer.getAtomPos(bond1->getBeginAtomIdx());
        auto end1_pos = conformer.getAtomPos(bond1->getEndAtomIdx());

        for (size_t j = i + 1; j < monomer_mol.getNumBonds(); j++) {
            auto bond2 = monomer_mol.getBondWithIdx(j);
            // skip if the bonds share an atom
            if (bond1->getBeginAtomIdx() == bond2->getBeginAtomIdx() ||
                bond1->getBeginAtomIdx() == bond2->getEndAtomIdx() ||
                bond1->getEndAtomIdx() == bond2->getBeginAtomIdx() ||
                bond1->getEndAtomIdx() == bond2->getEndAtomIdx()) {
                continue;
            }
            auto begin2_pos = conformer.getAtomPos(bond2->getBeginAtomIdx());
            auto end2_pos = conformer.getAtomPos(bond2->getEndAtomIdx());
            if (bonds_are_too_close(begin1_pos, end1_pos, begin2_pos,
                                    end2_pos)) {
                return false;
            }
        }
    }
    return true;
}

static bool has_no_stretched_bonds(const RDKit::ROMol& monomer_mol)
{
    const auto& conformer = monomer_mol.getConformer();
    for (auto bond : monomer_mol.bonds()) {
        float flexibility = 1.0;
        // non-backbone bonds can be more flexible so we allow them to stretch
        // more without flagging a failure
        if (!is_backbone_bond(bond)) {
            flexibility = NON_BACKBONE_BOND_FLEXIBILITY;
        }
        auto begin_pos = conformer.getAtomPos(bond->getBeginAtomIdx());
        auto end_pos = conformer.getAtomPos(bond->getEndAtomIdx());
        auto bond_length = (end_pos - begin_pos).length();
        if (bond_length >
            MAX_BOND_STRETCH * MONOMER_BOND_LENGTH * flexibility) {
            return false;
        }
    }
    return true;
}

bool coordinates_are_valid(const RDKit::ROMol& monomer_mol)
{
    return has_no_clashes(monomer_mol) && has_no_bond_crossings(monomer_mol) &&
           has_no_stretched_bonds(monomer_mol);
}

/**
 * Internal representations of a resize operation.
 */
struct MonomerResizeData {
    unsigned int index;         // Monomer index being resized
    RDGeom::Point3D ref_pos;    // Reference monomer position
    RDGeom::Point3D delta;      // Size delta (new - old)
    unsigned int x_cluster = 0; // Vertical stack identifier
    unsigned int y_cluster = 0; // Horizontal stack identifier
};

struct RingResizeData {
    unsigned int ring_id; // Ring index being resized
    double scale_factor;  // Coordinates scale factor (new / old)
};

/**
 * Assign cluster IDs along a given axis.
 *
 * @param resize_data  vector of monomers to cluster
 * @param coord_member pointer to either the x or y coordinate
 * @param cluster_member pointer to either the x_cluster or y_cluster member
 * @param eps clustering epsilon (distance threshold)
 */
template <typename TCoord, typename TCluster>
void assign_clusters(std::vector<MonomerResizeData>& resize_data,
                     TCoord RDGeom::Point3D::*coord_member,
                     TCluster MonomerResizeData::*cluster_member, double eps)
{
    unsigned int next_cluster_id = 1;

    for (size_t i = 0; i < resize_data.size(); ++i) {
        if (resize_data[i].*cluster_member != 0) {
            continue; // already assigned
        }

        resize_data[i].*cluster_member = next_cluster_id;

        for (size_t j = i + 1; j < resize_data.size(); ++j) {
            if (std::abs(resize_data[i].delta.*coord_member -
                         resize_data[j].delta.*coord_member) < eps) {
                resize_data[j].*cluster_member = next_cluster_id;
            }
        }

        ++next_cluster_id;
    }
}

/**
 * Read the stored monomer size for an atom.
 * Falls back to MONOMER_MINIMUM_SIZE when not present.
 */
inline RDGeom::Point3D get_monomer_size(const RDKit::ROMol& mol,
                                        unsigned int index)
{
    RDGeom::Point3D size(MONOMER_MINIMUM_SIZE, MONOMER_MINIMUM_SIZE, 0);

    mol.getAtomWithIdx(index)->getPropIfPresent<RDGeom::Point3D>(
        MONOMER_ITEM_SIZE, size);

    return size;
}

struct RingResizeInfo {
    std::vector<std::vector<int>> rings; // list of resizeable rings
    std::unordered_map<int, std::set<int>>
        atom_to_ring_map; // atom index -> set of ring indices
};

/**
 * Convert user resize requests for monomers in rings into internal resize
 * data. Compute how much each monomer's resize would contribute to the
 * expanding or shrinking of the ring, then store the maximum change as a scale
 * factor to apply to all the ring's monomer coordinates. This function performs
 * NO geometry mutation.
 */
std::vector<RingResizeData> collect_ring_resize_data(
    const RDKit::ROMol& mol,
    const std::unordered_map<int, RDGeom::Point3D>& resizes,
    const RingResizeInfo& ring_resize_info)
{
    std::vector<RingResizeData> result;
    auto& conformer = mol.getConformer();

    for (size_t ring_id = 0; ring_id < ring_resize_info.rings.size();
         ++ring_id) {
        const auto& ring = ring_resize_info.rings[ring_id];

        // --- centroid ---
        auto centroid = compute_centroid(mol, ring);

        // --- compute average current side length ---
        double side_sum = 0.0;
        for (size_t i = 0; i < ring.size(); ++i) {
            const auto& p0 = conformer.getAtomPos(ring[i]);
            const auto& p1 = conformer.getAtomPos(ring[(i + 1) % ring.size()]);
            side_sum += (p1 - p0).length();
        }
        const double current_side = side_sum / ring.size();
        if (current_side < EPSILON) {
            continue;
        }

        // --- find max projected side delta ---
        double max_side_delta = 0.0;

        for (const auto& atom_idx : ring) {
            auto resize_it = resizes.find(atom_idx);
            if (resize_it == resizes.end()) {
                continue;
            }

            const RDGeom::Point3D& new_size = resize_it->second;
            const RDGeom::Point3D current_size =
                get_monomer_size(mol, atom_idx);

            RDGeom::Point3D delta_size = new_size - current_size;

            RDGeom::Point3D pos = conformer.getAtomPos(atom_idx);
            RDGeom::Point3D radial = pos - centroid;

            radial.normalize();

            RDGeom::Point3D tangent(-radial.y, radial.x, 0.0);

            const double side_delta =
                delta_size.x * tangent.x + delta_size.y * tangent.y;

            max_side_delta = std::max(max_side_delta, side_delta);
        }

        if (max_side_delta > 0.0) {
            const double scale = (current_side + max_side_delta) / current_side;

            RingResizeData r;
            r.ring_id = ring_id;
            r.scale_factor = scale;
            result.push_back(r);
        }
    }

    return result;
}

/**
 * Convert user resize requests for monomers in chains into internal resize
 * data. This function performs NO geometry mutation.
 */
std::vector<MonomerResizeData> collect_linear_resize_data(
    const RDKit::ROMol& mol,
    const std::unordered_map<int, RDGeom::Point3D>& resizes,
    const RingResizeInfo& ring_resize_info)
{
    std::vector<MonomerResizeData> result;
    auto& conformer = mol.getConformer();

    result.reserve(resizes.size());

    for (const auto& r : resizes) {
        const unsigned int index = r.first;
        // if the monomer is in a ring, skip it
        if (ring_resize_info.atom_to_ring_map.contains(index)) {
            continue;
        }
        const RDGeom::Point3D& new_size = r.second;

        const RDGeom::Point3D current_size = get_monomer_size(mol, index);

        RDGeom::Point3D delta = new_size - current_size;
        // Skip no-op resizes
        if (delta.x == 0. && delta.y == 0.) {
            continue;
        }

        result.push_back({index, conformer.getAtomPos(index), delta});
    }

    // group vertically stacked monomers to avoid compound horizontal shifts
    assign_clusters(result, &RDGeom::Point3D::x, &MonomerResizeData::x_cluster,
                    MONOMER_MINIMUM_SIZE * 0.25);

    // group horizontally stacked monomers to avoid compound vertical shifts
    assign_clusters(result, &RDGeom::Point3D::y, &MonomerResizeData::y_cluster,
                    MONOMER_MINIMUM_SIZE * 0.25);

    return result;
}

/**
 * Compute per-monomer displacement vectors based on resize data.
 *
 * Horizontal and Vertical displacements are computed separately and clustered:
 *  - applied at most once per cluster
 *  - magnitude is the MAX delta among monomers in that cluster
 *  Rings are displaced as a rigid body from their centroid.

 * No atom positions are modified here.
 */
std::vector<RDGeom::Point3D>
compute_linear_displacements(const RDKit::ROMol& mol,
                             const std::vector<MonomerResizeData>& resize_data,
                             const RingResizeInfo& ring_resize_info)
{
    auto& conformer = mol.getConformer();
    const unsigned int num_atoms = conformer.getNumAtoms();

    std::vector<RDGeom::Point3D> displacements(num_atoms,
                                               RDGeom::Point3D(0, 0, 0));

    if (resize_data.empty()) {
        return displacements;
    }

    /**
     * Precompute maximum deltas per cluster for both axes
     */
    std::unordered_map<unsigned int, double> max_dx_per_cluster;
    std::unordered_map<unsigned int, double> max_dy_per_cluster;

    for (const auto& r : resize_data) {
        max_dx_per_cluster[r.x_cluster] =
            std::max(max_dx_per_cluster[r.x_cluster], r.delta.x);
        max_dy_per_cluster[r.y_cluster] =
            std::max(max_dy_per_cluster[r.y_cluster], r.delta.y);
    }

    /**
     * Compute displacement along a single axis for one atom,
     * applying each cluster at most once.
     */
    auto compute_axis_displacement =
        [](double pos, double ref, unsigned int cluster,
           const std::unordered_map<unsigned int, double>& max_per_cluster,
           std::unordered_set<unsigned int>& used_clusters) -> double {
        if (used_clusters.contains(cluster)) {
            return 0.0;
        }

        const double max_delta = max_per_cluster.at(cluster);

        if (pos > ref + MONOMER_MINIMUM_SIZE) {
            used_clusters.insert(cluster);
            return max_delta / 2.0;
        }
        if (pos < ref - MONOMER_MINIMUM_SIZE) {
            used_clusters.insert(cluster);
            return -max_delta / 2.0;
        }
        return 0.0;
    };

    /**
     * Compute displacement for a given position against all resize refs. For
     * monomers in chains this is used directly. For monomers in rings,
     * displacement is computed at the ring centroid and applied uniformly to
     * all ring monomers.
     */
    auto compute_displacement_for_position =
        [&](const RDGeom::Point3D& pos) -> RDGeom::Point3D {
        RDGeom::Point3D disp(0, 0, 0);

        std::unordered_set<unsigned int> used_x_clusters;
        std::unordered_set<unsigned int> used_y_clusters;

        for (const auto& r : resize_data) {
            disp.x +=
                compute_axis_displacement(pos.x, r.ref_pos.x, r.x_cluster,
                                          max_dx_per_cluster, used_x_clusters);

            disp.y +=
                compute_axis_displacement(pos.y, r.ref_pos.y, r.y_cluster,
                                          max_dy_per_cluster, used_y_clusters);
        }
        return disp;
    };

    /**
     * Linear (non-ring) atoms
     */
    for (unsigned int i = 0; i < num_atoms; ++i) {
        if (ring_resize_info.atom_to_ring_map.contains(i)) {
            continue;
        }
        const RDGeom::Point3D pos = conformer.getAtomPos(i);
        displacements[i] = compute_displacement_for_position(pos);
    }

    /**
     * Rings: compute displacement at centroid and apply uniformly
     */
    for (const auto& ring : ring_resize_info.rings) {
        auto centroid = compute_centroid(mol, ring);
        const RDGeom::Point3D ring_disp =
            compute_displacement_for_position(centroid);

        for (auto atom_idx : ring) {
            displacements[atom_idx] = ring_disp;
        }
    }
    return displacements;
}

std::set<int>
find_all_connected_monomers_outside_ring(const RDKit::ROMol& mol, int start_idx,
                                         const std::vector<int>& ring)
{
    std::set<int> result;
    std::set<int> visited;
    std::queue<int> to_visit;

    to_visit.push(start_idx);
    for (auto monomer_idx : ring) {
        visited.insert(monomer_idx);
    }
    visited.insert(start_idx);

    std::set<int> ring_set(ring.begin(), ring.end());

    while (!to_visit.empty()) {
        int current_idx = to_visit.front();
        to_visit.pop();

        for (auto neighbor :
             mol.atomNeighbors(mol.getAtomWithIdx(current_idx))) {
            int neighbor_idx = neighbor->getIdx();
            if (visited.contains(neighbor_idx)) {
                continue;
            }
            result.insert(neighbor_idx);
            visited.insert(neighbor_idx);
            to_visit.push(neighbor_idx);
        }
    }
    return result;
}

/**
 * Compute per-monomer displacement vectors based on ring expansion data.
 *
 * Rings are expanded radially outwards from their centroid, so that the regular
 * polygon shape is preserved.
 *
 * No atom positions are modified here.
 */
std::vector<RDGeom::Point3D> compute_ring_expansion_displacements(
    const RDKit::ROMol& mol, const std::vector<RingResizeData>& resize_data,
    const RingResizeInfo& ring_resize_info)
{
    auto& conformer = mol.getConformer();
    const unsigned int num_atoms = conformer.getNumAtoms();

    std::vector<RDGeom::Point3D> displacements(num_atoms,
                                               RDGeom::Point3D(0, 0, 0));

    if (resize_data.empty()) {
        return displacements;
    }

    for (const auto& r : resize_data) {
        const int ring_id = r.ring_id;
        const double scale = r.scale_factor;

        if (ring_id < 0 ||
            static_cast<size_t>(ring_id) >= ring_resize_info.rings.size()) {
            continue;
        }

        const auto& ring = ring_resize_info.rings[ring_id];

        auto centroid = compute_centroid(mol, ring);

        // --- scale ring ---
        for (auto monomer_idx : ring) {
            const RDGeom::Point3D pos = conformer.getAtomPos(monomer_idx);
            const RDGeom::Point3D new_pos =
                centroid + ((pos - centroid) * scale);
            auto monomer_displacement = new_pos - pos;
            // the displacement also applies to all monomers outside of the ring
            // that are connected to this monomer
            auto all_connected_monomers_outside_ring =
                find_all_connected_monomers_outside_ring(mol, monomer_idx,
                                                         ring);

            displacements[monomer_idx] += monomer_displacement;
            for (auto bound_monomer_idx : all_connected_monomers_outside_ring) {
                displacements[bound_monomer_idx] += monomer_displacement;
            }
        }
    }

    return displacements;
}

/**
 * Apply precomputed displacements to atom positions.
 * Each atom is moved at most once.
 */
void apply_displacements(RDKit::ROMol& mol,
                         const std::vector<RDGeom::Point3D>& displacements)
{
    auto& conformer = mol.getConformer();

    for (unsigned int i = 0; i < displacements.size(); ++i) {
        if (displacements[i].x == 0. && displacements[i].y == 0.) {
            continue;
        }
        auto pos = conformer.getAtomPos(i);
        pos += displacements[i];
        conformer.setAtomPos(i, pos);
    }
}

/**
 * Persist new monomer sizes on atoms.
 */
void update_monomer_sizes(
    RDKit::ROMol& mol,
    const std::unordered_map<int, RDGeom::Point3D>& monomer_sizes)
{
    for (const auto& r : monomer_sizes) {
        mol.getAtomWithIdx(r.first)->setProp<RDGeom::Point3D>(MONOMER_ITEM_SIZE,
                                                              r.second);
    }
}

bool is_geometrically_regular_ring_2d(const RDKit::ROMol& mol,
                                      const std::vector<int>& ring_atoms,
                                      double radius_tol, double angle_tol)
{
    const auto& conf = mol.getConformer();
    const size_t N = ring_atoms.size();
    if (N < 3) {
        return false;
    }

    auto centroid = compute_centroid(mol, ring_atoms);

    // --- radii + angles
    std::vector<double> radii;
    std::vector<double> angles;
    radii.reserve(N);
    angles.reserve(N);

    for (auto idx : ring_atoms) {
        const auto& p = conf.getAtomPos(idx);
        double dx = p.x - centroid.x;
        double dy = p.y - centroid.y;

        radii.push_back(std::hypot(dx, dy));
        angles.push_back(std::atan2(dy, dx));
    }

    // --- radius consistency
    double mean_r = std::accumulate(radii.begin(), radii.end(), 0.0) / N;

    double var = 0.0;
    for (double r : radii) {
        var += (r - mean_r) * (r - mean_r);
    }

    double stddev = std::sqrt(var / N);
    if (mean_r < EPSILON || stddev / mean_r > radius_tol) {
        return false;
    }

    // --- angular uniformity
    std::sort(angles.begin(), angles.end());

    const double expected = 2.0 * M_PI / N;
    double max_dev = 0.0;

    for (size_t i = 0; i < N; ++i) {
        double a0 = angles[i];
        double a1 = angles[(i + 1) % N];

        double delta = a1 - a0;
        if (i + 1 == N) {
            delta += 2.0 * M_PI; // wrap
        }

        max_dev = std::max(max_dev, std::abs(delta - expected));
    }

    return max_dev < angle_tol;
}

RingResizeInfo compute_ring_info_for_resize(RDKit::ROMol& mol)
{
    RingResizeInfo result;
    const auto& all_rings = mol.getRingInfo()->atomRings();

    std::vector<std::vector<int>> candidate_rings;

    // Step 1: filter by geometric regularity
    for (const auto& ring : all_rings) {
        if (!is_geometrically_regular_ring_2d(mol, ring)) {
            continue;
        }
        candidate_rings.push_back(ring);
    }
    for (const auto& ring : candidate_rings) {
        int ring_idx = static_cast<int>(result.rings.size());
        result.rings.push_back(ring);
        for (int atom_idx : ring) {
            result.atom_to_ring_map[atom_idx].insert(ring_idx);
        }
    }
    return result;
}

/**
 * Resize multiple monomers in a single, batched operation. The
 * operations differ for monomers in rings and in chains:
 *   - For monomers in rings, the ring is expanded radially outwards to
 * maintain its regular polygon shape.
 *   - For monomers in chains, the resize is batched to try and keep
 * vertical and horizontal alignments, minimizing compound shifts.
 *
 * @param mol           The molecule containing the monomers to resize.
 * @param monomer_sizes A map from atom indices to their new desired sizes.
 */
void resize_monomers(RDKit::ROMol& mol,
                     std::unordered_map<int, RDGeom::Point3D> monomer_sizes)
{
    if (monomer_sizes.empty()) {
        return;
    }
    // override ring info to ensure rings are fully perceived.
    compute_full_ring_info(mol);

    auto ring_resize_info = compute_ring_info_for_resize(mol);

    auto linear_resize_data =
        collect_linear_resize_data(mol, monomer_sizes, ring_resize_info);
    auto ring_resize_data =
        collect_ring_resize_data(mol, monomer_sizes, ring_resize_info);
    if (linear_resize_data.empty() && ring_resize_data.empty()) {
        return;
    }

    // compute displacements caused by monomers in chains
    auto linear_displacements =
        compute_linear_displacements(mol, linear_resize_data, ring_resize_info);

    // compute displacements caused by monomers in rings
    auto ring_displacements = compute_ring_expansion_displacements(
        mol, ring_resize_data, ring_resize_info);

    // apply movement
    apply_displacements(mol, linear_displacements);
    apply_displacements(mol, ring_displacements);

    // persist new sizes
    update_monomer_sizes(mol, monomer_sizes);

    // reset ring info to avoid leaking "internal" perception of rings
    mol.getRingInfo()->reset();
}

unsigned int compute_monomer_mol_coords(RDKit::ROMol& monomer_mol)
{
    // clear layout related props so we can start a fresh layout
    clear_layout_props(monomer_mol);
    auto [polymers, parent_polymer] = break_into_polymers(monomer_mol);
    // SHARED-9795: Special case for single polymers that are strictly
    // chains
    if (is_single_linear_polymer(monomer_mol)) {
        lay_out_snaked_linear_polymer(*polymers[0]);
    } else {
        lay_out_polymers(polymers, parent_polymer);
    }

    remove_cxsmiles_labels(monomer_mol);

    auto conformer_id =
        copy_polymer_coords_to_monomer_mol(monomer_mol, polymers);
    adjust_chem_polymer_coords(monomer_mol);

    // clear layout related props to prevent leaking "internal" props
    clear_layout_props(monomer_mol);

    if (!coordinates_are_valid(monomer_mol)) {
        // the coordinates are not good, try breaking the molecule into
        // topological units. This considers rings as a single unit, even if
        // they are made of monomers that belong to different polymers.
        // Branching chains are also considered separate units. create a copy of
        // monomer_mol to avoid modifying the original
        RDKit::ROMol monomer_mol_copy(monomer_mol);
        monomer_mol_copy.getRingInfo()->reset();
        auto [units, parent_unit] =
            break_into_topological_units(monomer_mol_copy);
        lay_out_polymers(units, parent_unit);
        monomer_mol_copy.clearConformers();
        copy_polymer_coords_to_monomer_mol(monomer_mol_copy, units);
        if (coordinates_are_valid(monomer_mol_copy)) {
            // the new coordinates are valid, copy them back to the original
            // mol. remove the previous conformer first
            monomer_mol.removeConformer(conformer_id);
            conformer_id =
                copy_polymer_coords_to_monomer_mol(monomer_mol, units);
        } else {
            std::cerr << "Warning: Generated coordinates for monomer mol "
                         "contain clashes, bond crossings, or stretched bonds."
                      << std::endl;
        }
        clear_layout_props(monomer_mol_copy);
    }

    return conformer_id;
}

} // namespace rdkit_extensions
} // namespace schrodinger
