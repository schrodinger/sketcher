// -------------------------------------------------------------------------
// Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/ChemTransforms/MolFragmenter.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>

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
constexpr double DIST_BETWEEN_MULTIPLE_POLYMERS = 5;
constexpr unsigned int MONOMERS_PER_SNAKE = 10;
const double PI = boost::math::constants::pi<double>();

constexpr double MONOMER_CLASH_DISTANCE = MONOMER_BOND_LENGTH * 0.25;
constexpr double BOND_CLASH_DISTANCE = MONOMER_CLASH_DISTANCE;
constexpr double MIN_ANGLE_BETWEEN_ADJACENT_BONDS =
    boost::math::constants::pi<double>() / 18; // 10 degrees

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
                // This is a connection attachment point
                polymer
                    .getBondBetweenAtoms(neighbor->getIdx(),
                                         monomer_to_place_idx)
                    ->hasProp(CUSTOM_BOND)) {
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

/**
 * Initializes RingInfo for `polymer` with all possible rings including any
 * rings that may have been formed with zero order bonds. We do this by
 * converting all the zero order bonds into single order bonds and then looking
 * for rings. (The bond types are restored to their original values by the end
 * of this method)
 */
static void compute_full_ring_info(const RDKit::ROMol& polymer)
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
 * Finds the largest ring in `polymer` that includes exactly one non-standard
 * connection (i.e. exactly one connection that's not a standard backbone
 * connection). This serves to exclude cycles that span multiple chains, such as
 * DNA base pairing. If no such ring exists, then return the largest ring
 * regardless of the number of non-standard connections.
 */
static std::vector<int> find_largest_ring(const RDKit::ROMol& polymer)
{
    if (!polymer.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(polymer, /*res=*/nullptr, include_dative_bonds);
    }

    std::vector<int> largest_ring_no_connections{};
    std::vector<int> largest_ring{};
    for (auto cycle : polymer.getRingInfo()->atomRings()) {
        if (cycle.size() <= largest_ring_no_connections.size()) {
            continue;
        }
        std::vector<unsigned int> connection_points{};
        for (size_t i = 0; i < cycle.size(); i++) {
            auto monomer = cycle[i];
            auto next_monomer = i == cycle.size() - 1 ? cycle[0] : cycle[i + 1];
            auto bond = polymer.getBondBetweenAtoms(monomer, next_monomer);

            if (bond->hasProp(CUSTOM_BOND)) {
                connection_points.push_back(i);
                if (connection_points.size() > 1) {
                    break;
                }
            }
        }
        if (connection_points.size() == 1) {
            // We want to rotate the ring vector so that the first monomer in
            // the connection bond is in the front. That way, the chains will
            // lay out directly to the right of the ring.
            largest_ring_no_connections = cycle;
            std::rotate(largest_ring_no_connections.begin(),
                        largest_ring_no_connections.begin() +
                            connection_points[0],
                        largest_ring_no_connections.end());
        } else if (cycle.size() > largest_ring.size()) {
            largest_ring = cycle;
            std::rotate(largest_ring.begin(),
                        largest_ring.begin() + connection_points[0],
                        largest_ring.end());
        }
    }
    if (largest_ring_no_connections.size() > 0) {
        return largest_ring_no_connections;
    }
    if (largest_ring.size() > 0) {
        return largest_ring;
    }

    throw std::runtime_error(
        "No cycles found with exactly one connection bond");
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

static void
lay_out_cyclic_polymer(RDKit::ROMol& polymer,
                       std::unordered_set<int>& placed_monomers_idcs)
{
    std::vector<int> largest_cycle = find_largest_ring(polymer);
    auto cycle_coords = get_coords_for_ngon(largest_cycle.size());
    auto& conformer = polymer.getConformer();

    // lay out all the monomers in the ring
    for (size_t i = 0; i < largest_cycle.size(); i++) {
        auto monomer = polymer.getAtomWithIdx(largest_cycle[i]);
        place_monomer_at(conformer, monomer, cycle_coords[i],
                         placed_monomers_idcs);
    }

    // lay out all the monomers attached to monomers in the ring
    for (auto monomer_idx : largest_cycle) {
        auto monomer = polymer.getAtomWithIdx(monomer_idx);
        for (auto neighbor : polymer.atomNeighbors(monomer)) {
            if (placed_monomers_idcs.contains(
                    neighbor->getProp<int>(ORIGINAL_INDEX))) {
                continue;
            }

            auto pos = conformer.getAtomPos(monomer_idx);
            // chains attached to monomers in the first bond of the cycle (which
            // is parallel to the y-axis) will just be laid out horizontally to
            // the right of the cycle
            if (monomer_idx == largest_cycle[0] ||
                monomer_idx == largest_cycle[1]) {
                pos.x = pos.x + MONOMER_BOND_LENGTH;
                // The other chains start at 1 MONOMER_BOND_LENGTH away along
                // the radius and then lay out horizontally from there
            } else {
                auto radius = pos.length();
                auto multiplier = (radius + MONOMER_BOND_LENGTH) / radius;
                pos = pos * multiplier;
            }

            // branch monomer i.e. single monomer branched off the polymer
            // backbone is just placed along the radius
            if (neighbor->getProp<bool>(BRANCH_MONOMER)) {
                place_monomer_at(conformer, neighbor, pos,
                                 placed_monomers_idcs);
            } else {
                ChainDirection chain_dir =
                    pos.x < 0 ? ChainDirection::RTL : ChainDirection::LTR;
                lay_out_chain(polymer, neighbor, placed_monomers_idcs, pos,
                              chain_dir);
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
 * Lays out a simple long linear polymer (with no branches) in a snaking pattern
 */
static void lay_out_snaked_linear_polymer(RDKit::ROMol& polymer)
{
    auto& conformer = polymer.getConformer();
    RDGeom::Point3D chain_start_pos(0, 0, 0);
    ChainDirection chain_dir = ChainDirection::LTR;
    auto placed_monomers_idcs = std::unordered_set<int>{};

    for (size_t i = 0; i < polymer.getNumAtoms(); i += MONOMERS_PER_SNAKE) {
        auto start_idx = i;
        auto end_idx = start_idx + MONOMERS_PER_SNAKE;

        lay_out_chain(polymer, polymer.getAtomWithIdx(start_idx),
                      placed_monomers_idcs, chain_start_pos, chain_dir,
                      BranchDirection::UP,
                      [end_idx](const RDKit::Atom* monomer) {
                          return monomer->getIdx() == end_idx;
                      });

        if (end_idx < polymer.getNumAtoms()) {
            // next chain will be under this one and in the opposite direction
            chain_dir = chain_dir == ChainDirection::LTR ? ChainDirection::RTL
                                                         : ChainDirection::LTR;
            chain_start_pos = conformer.getAtomPos(end_idx - 1);
            chain_start_pos.y -= MONOMER_BOND_LENGTH;
        }
    }
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
// Assign unit IDs to rings. takes a map of atom idx to unit ID string to fill
// it. Each atom in a ring gets the same unit ID.
static void assign_ring_unit_ids(const RDKit::ROMol& mol,
                                 std::unordered_map<int, std::string>& unit_ids)
{
    if (!mol.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(mol, /*res=*/nullptr, include_dative_bonds);
    }

    for (auto cycle : mol.getRingInfo()->atomRings()) {
        std::string ring_id = "ring_" + std::to_string(cycle[0]);
        for (auto atom_idx : cycle) {
            unit_ids[atom_idx] = ring_id;
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

                    // Assign chain IDs, longest gets same as ring, branches get
                    // new IDs

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
 * Translate `polymer` so that it is positioned below `reference`
 * without overlap, while preserving its internal geometry.
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

    // Vertical separation needed to avoid overlap with the reference polymer
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
 * to avoid geometric clashes between inter-polymer bonds (polymer <->
 * reference) Special handling is applied for adjacent bonds sharing an atom,
 * since they represent a chemically valid continuation and require
 * a more local proximity check.
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
        auto begin_pos = conformer.getAtomPos(bond->getBeginAtomIdx());
        auto end_pos = conformer.getAtomPos(bond->getEndAtomIdx());
        auto bond_length = (end_pos - begin_pos).length();
        if (bond_length > MAX_BOND_STRETCH * MONOMER_BOND_LENGTH) {
            return false;
        }
    }
    return true;
}

static bool coordinates_are_valid(RDKit::ROMol& monomer_mol)
{
    return has_no_clashes(monomer_mol) && has_no_bond_crossings(monomer_mol) &&
           has_no_stretched_bonds(monomer_mol);
}

void resize_monomer(RDKit::ROMol& monomer_mol, unsigned int index,
                    const RDGeom::Point3D& new_size)
{
    // get current size
    auto atom = monomer_mol.getAtomWithIdx(index);
    RDGeom::Point3D current_size(MONOMER_MINIMUM_SIZE, MONOMER_MINIMUM_SIZE, 0);
    atom->getPropIfPresent<RDGeom::Point3D>(MONOMER_ITEM_SIZE, current_size);

    // calculate difference
    auto difference = new_size - current_size;
    if (difference.x == 0. && difference.y == 0.) {
        return;
    }

    // move every atom position accordingly
    auto& conformer = monomer_mol.getConformer();
    auto reference_monomer_position = conformer.getAtomPos(index);
    for (unsigned int i = 0u; i < conformer.getNumAtoms(); ++i) {
        if (i == index) {
            continue;
        }
        auto atom_pos = conformer.getAtomPos(i);
        if (atom_pos.x > reference_monomer_position.x + MONOMER_MINIMUM_SIZE) {
            atom_pos.x += difference.x / 2;
        } else if (atom_pos.x <
                   reference_monomer_position.x - MONOMER_MINIMUM_SIZE) {
            atom_pos.x -= difference.x / 2;
        }
        if (atom_pos.y > reference_monomer_position.y + MONOMER_MINIMUM_SIZE) {
            atom_pos.y += difference.y / 2;
        } else if (atom_pos.y <
                   reference_monomer_position.y - MONOMER_MINIMUM_SIZE) {
            atom_pos.y -= difference.y / 2;
        }
        conformer.setAtomPos(i, atom_pos);
    }
    atom->setProp<RDGeom::Point3D>(MONOMER_ITEM_SIZE, new_size);
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
