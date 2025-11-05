// -------------------------------------------------------------------------
// Copyright Schrodinger LLC, All Rights Reserved.
#include "schrodinger/rdkit_extensions/helm.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/ChemTransforms/MolFragmenter.h>
#include <rdkit/Geometry/point.h>
#include <boost/geometry.hpp>

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
constexpr double MONOMER_BOND_LENGTH = 1.5;
constexpr double DIST_BETWEEN_MULTIPLE_POLYMERS = 5;
constexpr unsigned int MONOMERS_PER_SNAKE = 10;
const double PI = boost::math::constants::pi<double>();

constexpr double MONOMER_CLASH_DISTANCE = MONOMER_BOND_LENGTH * 0.25;
constexpr double BOND_CLASH_DISTANCE = MONOMER_CLASH_DISTANCE;

// Props used in the fragmented polymer mols to store monomer graph info from
// the monomersitc mol
const std::string BOND_TO{"bondTo"};
const std::string ORIGINAL_INDEX{"originalIndex"};
// Used to indicate whether a monomer has had its position calculated and set in
// the conformer
const std::string MONOMER_PLACED{"monomerPlaced"};
// Set on the polymer mol
const std::string POLYMER_ID{"polymerID"};

// Replacement atom label for monomers with SMILES strings as their atom label
const std::string SMILES_MONOMER_LABEL{"CX"};

enum ChainDirection { LTR, RTL };
enum BranchDirection { UP, DOWN };
enum PolygonStartSide { LEFT, RIGHT };

namespace bg = boost::geometry;
using point_t = bg::model::d2::point_xy<double>;
using segment_t = bg::model::segment<point_t>;

static auto default_stop_condition = [](const RDKit::Atom* monomer) {
    return false;
};
/**
 * Lays out a monomer chain in `polymer` by traversing through the monomer graph
 * until either the `stop_condition` is met for all visited monomers or we run
 * out of unplaced monomers reachable from start_monomer. Expects all monomers
 * to have the MONOMER_PLACED prop set to indicate if they've already been
 * placed (polymers returned from `break_into_polymers` below already have this
 * prop set on all monomers).
 */
template <typename _Cond = decltype(default_stop_condition)> static void
lay_out_chain(RDKit::ROMol& polymer, const RDKit::Atom* start_monomer,
              const RDGeom::Point3D& start_pos = RDGeom::Point3D(0, 0, 0),
              ChainDirection chain_dir = LTR,
              BranchDirection branch_direction = UP,
              _Cond stop_condition = default_stop_condition)
{
    auto& conformer = polymer.getConformer();
    auto x_pos = start_pos.x;
    auto monomer_to_layout = start_monomer;
    while (monomer_to_layout != nullptr) {
        if (stop_condition(monomer_to_layout)) {
            break;
        }
        auto monomer_idx = monomer_to_layout->getIdx();
        conformer.setAtomPos(monomer_idx,
                             RDGeom::Point3D(x_pos, start_pos.y, start_pos.z));
        monomer_to_layout->setProp(MONOMER_PLACED, true);

        auto monomer_neighbors = polymer.atomNeighbors(monomer_to_layout);
        monomer_to_layout = nullptr;
        for (auto neighbor : monomer_neighbors) {
            if (neighbor->getProp<bool>(MONOMER_PLACED) ||
                // This is a connection attachment point
                polymer.getBondBetweenAtoms(neighbor->getIdx(), monomer_idx)
                    ->hasProp(CUSTOM_BOND)) {
                continue;
            }

            if (!neighbor->getProp<bool>(BRANCH_MONOMER)) {
                monomer_to_layout = neighbor;
            } else {
                // branch monomer - placed with a y offset at the same x coord
                double y_offset = branch_direction == UP ? MONOMER_BOND_LENGTH
                                                         : -MONOMER_BOND_LENGTH;
                RDGeom::Point3D pos(x_pos, start_pos.y + y_offset, start_pos.z);
                conformer.setAtomPos(neighbor->getIdx(), pos);
                neighbor->setProp(MONOMER_PLACED, true);
            }
        }
        x_pos += chain_dir == LTR ? MONOMER_BOND_LENGTH : -MONOMER_BOND_LENGTH;
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
 * Finds the largest ring in `polymer` with exactly one connection attachment
 * point
 */
static std::vector<int> find_largest_ring(const RDKit::ROMol& polymer)
{
    if (!polymer.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(polymer, /*res=*/nullptr, include_dative_bonds);
    }

    std::vector<int> largest_ring{};
    for (auto cycle : polymer.getRingInfo()->atomRings()) {
        if (cycle.size() <= largest_ring.size()) {
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
            largest_ring = cycle;
            std::rotate(largest_ring.begin(),
                        largest_ring.begin() + connection_points[0],
                        largest_ring.end());
        }
    }

    if (largest_ring.size() == 0) {
        throw std::runtime_error(
            "No cycles found with exactly one connection bond");
    }

    return largest_ring;
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
get_coords_for_ngon(size_t n, PolygonStartSide start_side = RIGHT)
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

static void lay_out_cyclic_polymer(RDKit::ROMol& polymer)
{
    std::vector<int> largest_cycle = find_largest_ring(polymer);
    auto cycle_coords = get_coords_for_ngon(largest_cycle.size());
    auto& conformer = polymer.getConformer();

    // lay out all the monomers in the ring
    for (size_t i = 0; i < largest_cycle.size(); i++) {
        auto monomer = polymer.getAtomWithIdx(largest_cycle[i]);
        conformer.setAtomPos(monomer->getIdx(), cycle_coords[i]);
        monomer->setProp(MONOMER_PLACED, true);
    }

    // lay out all the monomers attached to monomers in the ring
    for (auto monomer_idx : largest_cycle) {
        auto monomer = polymer.getAtomWithIdx(monomer_idx);
        for (auto neighbor : polymer.atomNeighbors(monomer)) {
            if (neighbor->getProp<bool>(MONOMER_PLACED)) {
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
                conformer.setAtomPos(neighbor->getIdx(), pos);
                neighbor->setProp(MONOMER_PLACED, true);
            } else {
                ChainDirection chain_dir = pos.x < 0 ? RTL : LTR;
                lay_out_chain(polymer, neighbor, pos, chain_dir);
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
static void layout_hairpin_polymer(RDKit::ROMol& polymer)
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
        conformer.setAtomPos(monomer->getIdx(), coords);
        monomer->setProp(MONOMER_PLACED, true);

        for (auto neighbor : polymer.atomNeighbors(monomer)) {
            if (neighbor->getProp<bool>(BRANCH_MONOMER)) {
                auto radius = coords.length();
                auto multiplier = (radius + MONOMER_BOND_LENGTH) / radius;
                conformer.setAtomPos(neighbor->getIdx(), coords * multiplier);
                neighbor->setProp(MONOMER_PLACED, true);
            }
        }
    }

    // layout the chains
    lay_out_chain(polymer, hairpin_turn.front(), cycle_coords[turn_start],
                  ChainDirection::RTL, BranchDirection::DOWN);
    lay_out_chain(polymer, hairpin_turn.back(), cycle_coords[turn_end],
                  ChainDirection::RTL, BranchDirection::UP);
}

static void lay_out_linear_polymer(RDKit::ROMol& polymer,
                                   const bool rotate = false)
{
    auto monomers = polymer.atoms();
    auto start_monomer =
        *std::find_if(monomers.begin(), monomers.end(), [](RDKit::Atom* m) {
            return !m->getProp<bool>(BRANCH_MONOMER);
        });
    auto branch_dir = rotate ? BranchDirection::UP : BranchDirection::DOWN;
    auto chain_dir = rotate ? ChainDirection::RTL : ChainDirection::LTR;
    lay_out_chain(polymer, start_monomer, RDGeom::Point3D(0, 0, 0), chain_dir,
                  branch_dir);
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
static void lay_out_polymer(RDKit::ROMol& polymer, const bool rotate = false)
{
    if (!polymer.getRingInfo()->isInitialized()) {
        constexpr bool include_dative_bonds = true;
        RDKit::MolOps::findSSSR(polymer, /*res=*/nullptr, include_dative_bonds);
    }
    if (polymer.getRingInfo()->numRings() > 0) {
        lay_out_cyclic_polymer(polymer);
    } else if (polymer.getNumAtoms() == polymer.getNumBonds() + 1) {
        lay_out_linear_polymer(polymer, rotate);
    } else {
        layout_hairpin_polymer(polymer);
    }
}

/**
 * Lays out a simple long linear polymer (with no branches) in a snaking pattern
 */
static void lay_out_snaked_linear_polymer(RDKit::ROMol& polymer)
{
    auto& conformer = polymer.getConformer();
    RDGeom::Point3D chain_start_pos(0, 0, 0);
    ChainDirection chain_dir = LTR;

    for (size_t i = 0; i < polymer.getNumAtoms(); i += MONOMERS_PER_SNAKE) {
        auto start_idx = i;
        auto end_idx = start_idx + MONOMERS_PER_SNAKE;

        lay_out_chain(polymer, polymer.getAtomWithIdx(start_idx),
                      chain_start_pos, chain_dir, BranchDirection::UP,
                      [end_idx](const RDKit::Atom* monomer) {
                          return monomer->getIdx() == end_idx;
                      });

        if (end_idx < polymer.getNumAtoms()) {
            // next chain will be under this one and in the opposite direction
            chain_dir = chain_dir == LTR ? RTL : LTR;
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
        atom->setProp(MONOMER_PLACED, false);
    }

    for (auto bond : monomer_mol.bonds()) {
        auto beginMonomer = bond->getBeginAtom();
        auto endMonomer = bond->getEndAtom();
        if (get_polymer_id(beginMonomer) == get_polymer_id(endMonomer)) {
            continue;
        }
        // store the connected indices as comman separated values in the
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

/**
 * Provides the y-coordinate offset to apply to each monomer in
 * `polymer_to_translate` to place it one `MONOMER_BOND_LENGTH` away below
 * `reference_polymer` with no overlap
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
 * Move polymer_to_orient so that it is positioned below reference_polymer
 * with no overlap.
 */
static void orient_polymer(RDKit::ROMol& polymer_to_orient,
                           const RDKit::ROMol& reference_polymer)
{
    auto polymer_bonds =
        get_bonds_between_polymers(polymer_to_orient, reference_polymer);
    auto pos_offset = RDGeom::Point3D(0, 0, 0);
    if (polymer_bonds.empty()) {
        auto y_offset =
            get_y_offset_for_polymer(polymer_to_orient, reference_polymer);
        pos_offset.y = y_offset - DIST_BETWEEN_MULTIPLE_POLYMERS;
    } else {
        auto& monomer_coords =
            polymer_to_orient.getConformer().getAtomPos(polymer_bonds[0].first);
        auto& placed_monomer_coords =
            reference_polymer.getConformer().getAtomPos(
                polymer_bonds[0].second);
        auto y_offset =
            get_y_offset_for_polymer(polymer_to_orient, reference_polymer);
        // This is where the monomer on this polymer should end up after
        // translation - directly under the placed monomer so the bond
        // is vertical
        auto translated_monomer_coords = RDGeom::Point3D(
            placed_monomer_coords.x, monomer_coords.y + y_offset,
            placed_monomer_coords.z);
        pos_offset = translated_monomer_coords - monomer_coords;
    }
    for (auto& pos : polymer_to_orient.getConformer().getPositions()) {
        pos += pos_offset;
    }
}

/**
 * Adds a conformer to the monomer_mol copying all the coordinates from the
 * polymers. Expects each monomer in the polymers to have an ORIGINAL_INDEX
 * prop pointing to the corresponding index for that monomer in the
 * monomer_mol. Returns the id of the added conformer.
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
 * connect a polymer back to itself or two polymers together) to position
 * them between the two polymers they are connecting.
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
    std::vector<std::string> layout_props{BOND_TO, MONOMER_PLACED,
                                          ORIGINAL_INDEX};
    for (auto monomer : monomer_mol.atoms()) {
        for (auto prop : layout_props) {
            monomer->clearProp(prop);
        }
    }
}

/**
 * Determines if a monomer mol contains a single linear polymer with no
 * branches or rings
 */
static bool is_single_linear_polymer(const RDKit::ROMol& monomer_mol)
{
    auto polymer_ids = get_polymer_ids(monomer_mol);
    if (polymer_ids.size() != 1) {
        return false;
    }

    compute_full_ring_info(monomer_mol);
    if (monomer_mol.getRingInfo()->numRings() > 0) {
        // reset ring info so it can get recomputed by the appropriate
        // layout method
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
    // return true if the two polymers are nucleic acids and there are at
    // least two bonds between them.

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
    // lay out the polymers in connection order so connected polymers are
    // laid out next to each other.
    for (auto polymer : polymers) {
        RDKit::ROMOL_SPTR parent = nullptr;
        if (parent_polymer.contains(polymer)) {
            parent = parent_polymer.at(polymer);
        }
        // For double stranded nucleic acids we want to lay out the
        // first polymer normally and then rotate the other polymer
        // 180° so the strands run anti-parallel to each other
        bool rotate_polymer = (parent != nullptr) &&
                              are_double_stranded_nucleic_acid(polymer, parent);
        lay_out_polymer(*polymer, rotate_polymer);
        if (parent != nullptr) {
            orient_polymer(*polymer, *parent);
        }
    }
}

bool has_no_clashes(const RDKit::ROMol& monomer_mol)
{
    auto& conformer = monomer_mol.getConformer();
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

double compute_distance_between_segments(const RDGeom::Point3D& p1,
                                         const RDGeom::Point3D& p2,
                                         const RDGeom::Point3D& q1,
                                         const RDGeom::Point3D& q2)
{
    segment_t seg_p = {{p1.x, p1.y}, {p2.x, p2.y}};
    segment_t seg_q = {{q1.x, q1.y}, {q2.x, q2.y}};
    return bg::distance(seg_p, seg_q);
}

bool has_no_bond_crossings(const RDKit::ROMol& monomer_mol)
{
    auto& conformer = monomer_mol.getConformer();
    auto positions = conformer.getPositions();
    for (size_t i = 0; i < monomer_mol.getNumBonds(); i++) {
        auto bond1 = monomer_mol.getBondWithIdx(i);
        auto begin1_pos = conformer.getAtomPos(bond1->getBeginAtomIdx());
        auto end1_pos = conformer.getAtomPos(bond1->getEndAtomIdx());
        auto minx1 = std::min(begin1_pos.x, end1_pos.x) - BOND_CLASH_DISTANCE;
        auto maxx1 = std::max(begin1_pos.x, end1_pos.x) + BOND_CLASH_DISTANCE;
        auto miny1 = std::min(begin1_pos.y, end1_pos.y) - BOND_CLASH_DISTANCE;
        auto maxy1 = std::max(begin1_pos.y, end1_pos.y) + BOND_CLASH_DISTANCE;

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
            auto minx2 =
                std::min(begin2_pos.x, end2_pos.x) - BOND_CLASH_DISTANCE;
            auto maxx2 =
                std::max(begin2_pos.x, end2_pos.x) + BOND_CLASH_DISTANCE;
            auto miny2 =
                std::min(begin2_pos.y, end2_pos.y) - BOND_CLASH_DISTANCE;
            auto maxy2 =
                std::max(begin2_pos.y, end2_pos.y) + BOND_CLASH_DISTANCE;
            // bounding box check
            if (maxx1 < minx2 || maxx2 < minx1 || maxy1 < miny2 ||
                maxy2 < miny1) {
                continue;
            }
            if (segments_intersect(begin1_pos, end1_pos, begin2_pos,
                                   end2_pos) ||
                compute_distance_between_segments(begin1_pos, end1_pos,
                                                  begin2_pos, end2_pos) <
                    BOND_CLASH_DISTANCE) {
                return false;
            }
        }
    }
    return true;
}

[[maybe_unused]] static bool coordinates_are_valid(RDKit::ROMol& monomer_mol)
{
    return has_no_clashes(monomer_mol) && has_no_bond_crossings(monomer_mol);
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
    return conformer_id;
}

} // namespace rdkit_extensions
} // namespace schrodinger
