#include "schrodinger/sketcher/rdkit/fragment.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <tuple>

#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/QueryOps.h>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Get information about the fragment
 * @param fragment A molecule with exactly one attachment point
 * @return A tuple of
 *   - the fragment's attachment point parent atom, i.e., the heavy atom that
 *     the attachment point dummy atom is bound to
 *   - the fragment's attachment point dummy atom
 *   - whether the attachment point parent atom is in a ring
 */
std::tuple<const RDKit::Atom*, const RDKit::Atom*, bool>
get_frag_info(const RDKit::ROMol& fragment)
{
    std::vector<RDKit::Atom*> all_aps;
    auto all_atoms = fragment.atoms();
    std::copy_if(all_atoms.begin(), all_atoms.end(),
                 std::back_inserter(all_aps), is_attachment_point);
    if (all_aps.size() != 1) {
        throw std::runtime_error(
            "Fragments must have exactly one attachment point.");
    }
    auto* ap_dummy_atom = all_aps.front();
    auto* ap_parent_atom = *fragment.atomNeighbors(ap_dummy_atom).begin();
    bool is_ring = RDKit::queryIsAtomInRing(ap_parent_atom);
    return {ap_parent_atom, ap_dummy_atom, is_ring};
}

/**
 * Move the fragment so that it is in the appropriate position relative to the
 * core
 * @param fragment The fragment to overlay.  This molecule must have exactly one
 * attachment point.
 * @param frag_atom The fragment atom to overlay on top of core_atom
 * @param frag_bond The fragment bond to align with core_bond
 * @param core The molecule to overlay fragment on
 * @param core_atom The atom to overlay frag_atom on
 * @param core_bond The bond to align core_bond with
 * @return A new fragment conformer.  Note that this conformer is not added to
 * fragment.
 */
RDKit::Conformer overlay_fragment_on_core(const RDKit::ROMol& fragment,
                                          const RDKit::Atom* const frag_atom,
                                          const RDKit::Bond* const frag_bond,
                                          const RDKit::ROMol& core,
                                          const RDKit::Atom* const core_atom,
                                          const RDKit::Bond* const core_bond)
{
    const auto& core_conf = core.getConformer();
    const auto core_atom_idx = core_atom->getIdx();
    const auto& core_atom_coords = core_conf.getAtomPos(core_atom_idx);
    auto frag_conf = translate_fragment_ap_to(fragment, core_atom_coords);

    const auto frag_atom_idx = frag_atom->getIdx();
    const auto other_frag_atom_idx = frag_bond->getOtherAtomIdx(frag_atom_idx);
    const auto& other_frag_atom_coords =
        frag_conf.getAtomPos(other_frag_atom_idx);

    const auto other_core_atom_idx = core_bond->getOtherAtomIdx(core_atom_idx);
    const auto& other_core_atom_coords =
        core_conf.getAtomPos(other_core_atom_idx);
    auto rotate_by = get_angle_radians(other_frag_atom_coords, core_atom_coords,
                                       other_core_atom_coords);
    rotate_conformer_radians(rotate_by, core_atom_coords, frag_conf);
    return frag_conf;
}

/**
 * Determine whether a ring fragment should be flipped 180 degrees to avoid
 * clashing with core atoms.  To do this, we look at the core atoms that are
 * neighbors of the atoms involved in core_bond.  We want these neighbors to be
 * on the *same* side of core_bond as frag_ap_dummy_atom, since we assume that
 * frag_ap_dummy_atom is outside of the ring (which means that this would put
 * the ring on the *opposite* side of the bond from the neighbors.
 * @param fragment The fragment.  This molecule must have exactly one attachment
 * point and it must be bound to a ring atom.
 * @param frag_conf The conformer that overlays the fragment on the core
 * @param frag_ap_dummy_atom The fragment attachment point atom
 * @param core The core molecule
 * @param core_bond The core bond that the fragment is overlayed on top of
 * @return Whether the fragment should be flipped
 */
bool should_flip_ring_fragment(const RDKit::ROMol& fragment,
                               const RDKit::Conformer& frag_conf,
                               const RDKit::Atom* const frag_ap_dummy_atom,
                               const RDKit::ROMol& core,
                               const RDKit::Bond* const core_bond)
{
    unsigned int num_same_side = 0;
    unsigned int num_opposite_side = 0;
    auto& core_conf = core.getConformer();
    auto& frag_ap_coords = frag_conf.getAtomPos(frag_ap_dummy_atom->getIdx());
    for (const auto* atom :
         {core_bond->getBeginAtom(), core_bond->getEndAtom()}) {
        const auto* other_atom = core_bond->getOtherAtom(atom);
        const auto& atom_coords = core_conf.getAtomPos(atom->getIdx());
        const auto& other_atom_coords =
            core_conf.getAtomPos(other_atom->getIdx());
        auto relative_frag_ap_coords = frag_ap_coords - atom_coords;
        auto relative_other_atom_coords = other_atom_coords - atom_coords;

        for (const RDKit::Atom* neighbor : core.atomNeighbors(atom)) {
            if (neighbor == other_atom) {
                continue;
            }
            auto& neighbor_coords = core_conf.getAtomPos(neighbor->getIdx());
            auto relative_neighbor_coords = neighbor_coords - atom_coords;

            if (are_points_on_same_side_of_line(relative_neighbor_coords,
                                                relative_frag_ap_coords,
                                                relative_other_atom_coords)) {
                ++num_same_side;
            } else {
                ++num_opposite_side;
            }
        }
    }
    return num_opposite_side > num_same_side;
}

/**
 * Get a list of candidate bonds from the fragment that could possibly be used
 * for aligning the fragment with the core.  The bonds are all between the
 * attachment point parent atom and a non-attachment-point atom.  If the
 * parent atom is in a ring, then the bonds will all be ring bonds.
 * @param fragment The fragment molecule
 * @param ap_parent_atom The attachment point parent atom
 * @param is_ring Whether the attachment point parent atom is part of a ring
 * @return A list of all candidate bonds
 */
std::vector<const RDKit::Bond*>
get_candidate_frag_bonds_for_alignment(const RDKit::ROMol& fragment,
                                       const RDKit::Atom* const ap_parent_atom,
                                       const bool is_ring)
{
    // make a list of all of the fragment candidate bonds
    std::vector<const RDKit::Bond*> frag_bonds;
    auto frag_bonds_it = fragment.atomBonds(ap_parent_atom);

    std::function<bool(const RDKit::Bond* const)> bond_criterion;
    if (is_ring) {
        // if the fragment is a ring fragment, then we always want to use one of
        // the ring bonds
        bond_criterion = RDKit::queryIsBondInRing;
    } else {
        bond_criterion = [&ap_parent_atom](const RDKit::Bond* const bond) {
            auto other_atom = bond->getOtherAtom(ap_parent_atom);
            return !is_attachment_point(other_atom);
        };
    }
    std::copy_if(frag_bonds_it.begin(), frag_bonds_it.end(),
                 std::back_inserter(frag_bonds), bond_criterion);
    return frag_bonds;
}

/**
 * Determine which fragment bond should be aligned to the given core bond
 * @param fragment The fragment molecule
 * @param frag_ap_parent_atom The attachment point parent atom from the fragment
 * @param core The core molecule
 * @param core_bond The core bond to align to
 * @return The appropriate fragment bond
 */
const RDKit::Bond* get_frag_bond_to_align_to(
    const RDKit::ROMol& fragment, const RDKit::Atom* const frag_ap_parent_atom,
    const RDKit::ROMol& core, const RDKit::Bond* const core_bond)
{
    bool frag_is_ring = RDKit::queryIsAtomInRing(frag_ap_parent_atom);
    auto frag_bonds = get_candidate_frag_bonds_for_alignment(
        fragment, frag_ap_parent_atom, frag_is_ring);

    if (frag_is_ring && RDKit::queryIsBondInRing(core_bond)) {
        // if the fragment and the core are both rings, then frag_bond should be
        // the bond with the highest order
        return *std::max_element(
            frag_bonds.begin(), frag_bonds.end(),
            [](const RDKit::Bond* const bond1, const RDKit::Bond* const bond2) {
                return bond1->getBondTypeAsDouble() <
                       bond2->getBondTypeAsDouble();
            });
    } else {
        // Otherwise, we use the bond with the bond order (e.g. single vs.
        // double) closest to core_bond
        auto core_bond_order = core_bond->getBondTypeAsDouble();
        auto bond_degree_diff =
            [core_bond_order](const RDKit::Bond* const bond) {
                return std::abs(bond->getBondTypeAsDouble() - core_bond_order);
            };
        return *std::min_element(
            frag_bonds.begin(), frag_bonds.end(),
            [&bond_degree_diff](const RDKit::Bond* const bond1,
                                const RDKit::Bond* const bond2) {
                return bond_degree_diff(bond1) < bond_degree_diff(bond2);
            });
    }
}

/**
 * Return a non-dummy atom bound to an attachment point parent atom
 * @param mol The molecule
 * @param ap_parent_atom The attachment point parent atom (i.e. the heavy atom
 * bound to the attachment point dummy atom)
 * @param is_ring Whether ap_parent_atom is part of a ring.  If it is, then the
 * returned atom will always be a ring atom.
 */
const RDKit::Atom*
get_atom_bound_to_ap_parent_atom(const RDKit::ROMol& mol,
                                 const RDKit::Atom* const ap_parent_atom,
                                 const bool is_ring)
{
    auto bonds =
        get_candidate_frag_bonds_for_alignment(mol, ap_parent_atom, is_ring);
    auto* bond = bonds.front();
    return bond->getOtherAtom(ap_parent_atom);
}

/**
 * Return true if the atom is a heteroatom (not C or H) or a query.  Otherwise,
 * return false.
 */
bool is_heteroatom_or_query(const RDKit::Atom* const atom)
{
    if (atom->hasQuery()) {
        return true;
    }
    auto core_atomic_num = atom->getAtomicNum();
    return core_atomic_num != static_cast<int>(Element::C) &&
           core_atomic_num != static_cast<int>(Element::H);
}

} // namespace

static std::unordered_set<RDKit::Bond*>
determine_fragment_bonds_to_mutate_to_single_bonds(
    RDKit::ROMol& fragment, const RDKit::ROMol& core,
    const AtomToAtomMap& frag_atom_to_core_atom);

static std::tuple<const RDKit::Atom*, const RDKit::Atom*, const RDKit::Bond*>
get_frag_neighbor_and_corresponding_core_neighbor_and_bond(
    const RDKit::Atom* const frag_atom, const RDKit::Bond* const frag_bond,
    const RDKit::Atom* const core_atom,
    const AtomToAtomMap& frag_atom_to_core_atom);

static bool should_replace_core_bond(const RDKit::Bond* const frag_bond,
                                     const RDKit::Bond* const core_bond);

static AtomToAtomMap get_fragment_to_core_atom_map(
    const RDKit::ROMol& fragment, const RDKit::Conformer& frag_conf,
    const RDKit::Atom* frag_start_atom, const RDKit::ROMol& core,
    const RDKit::Atom* core_start_atom);

static AtomToAtomMap get_fragment_to_core_atom_map(
    const RDKit::ROMol& fragment, const RDKit::Atom* frag_start_atom,
    const RDKit::ROMol& core, const RDKit::Atom* core_start_atom);

static bool should_replace_core_atom(const RDKit::Atom* const frag_atom,
                                     const RDKit::Atom* const core_atom);

static std::vector<std::pair<unsigned int, AtomFunc>>
determine_core_atom_mutations(const AtomToAtomMap& frag_atom_to_core_atom);

static std::tuple<AtomPtrToFragBondMap,
                  std::vector<std::pair<unsigned int, BondFunc>>,
                  std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>>
determine_core_bond_additions_and_mutations(
    const RDKit::ROMol& core, const RDKit::ROMol& fragment,
    AtomToAtomMap frag_atom_to_core_atom);

static AtomIdxToFragBondMap convert_bond_map_from_ptrs_to_idxs(
    const RDKit::ROMol& core,
    const AtomPtrToFragBondMap& core_to_frag_bonds_by_ptr);

/**
 * Make a copy of the conformer and add the offset to it
 */
static RDKit::Conformer get_conformer_plus_offset(RDKit::Conformer conf,
                                                  const RDGeom::Point3D& offset)
{
    // note that this makes a copy of the conformer
    for (auto& coords : conf.getPositions()) {
        coords += offset;
    }
    return conf;
}

RDKit::Conformer translate_fragment_ap_to(const RDKit::ROMol& fragment,
                                          const RDGeom::Point3D& point)
{
    auto [frag_ap_parent_atom, frag_ap_dummy_atom, is_ring] =
        get_frag_info(fragment);
    auto& conf = fragment.getConformer();
    auto offset = point - conf.getAtomPos(frag_ap_parent_atom->getIdx());
    return get_conformer_plus_offset(conf, offset);
}

RDKit::Conformer translate_fragment_center_to(const RDKit::ROMol& fragment,
                                              const RDGeom::Point3D& point)
{
    const RDKit::Atom* frag_ap_parent_atom;
    const RDKit::Atom* frag_ap_dummy_atom;
    bool is_ring;
    std::tie(frag_ap_parent_atom, frag_ap_dummy_atom, is_ring) =
        get_frag_info(fragment);
    std::unordered_set<const RDKit::Atom*> non_ap_atoms;
    auto all_atoms_it = fragment.atoms();
    std::copy_if(all_atoms_it.begin(), all_atoms_it.end(),
                 std::inserter(non_ap_atoms, non_ap_atoms.begin()),
                 [frag_ap_dummy_atom](auto* cur_atom) {
                     return cur_atom != frag_ap_dummy_atom;
                 });
    auto offset = point - find_centroid(fragment, non_ap_atoms);
    return get_conformer_plus_offset(fragment.getConformer(), offset);
}

RDKit::Conformer align_fragment_with_atom(const RDKit::ROMol& fragment,
                                          const RDKit::Atom* const core_atom)
{
    const auto& core = core_atom->getOwningMol();
    const auto& core_conf = core.getConformer();
    const auto core_atom_idx = core_atom->getIdx();
    const auto& core_atom_coords = core_conf.getAtomPos(core_atom_idx);

    auto [frag_ap_parent_atom, frag_ap_dummy_atom, is_ring] =
        get_frag_info(fragment);

    auto core_neighbors = get_relative_positions_of_atom_neighbors(core_atom);
    auto core_offset = best_placing_around_origin(
        core_neighbors, /* limit_to_120_for_single_neighbor = */ !is_ring);
    auto core_best_placing = core_offset + core_atom_coords;

    const RDKit::Atom* frag_other_atom;
    if (is_ring) {
        frag_other_atom = frag_ap_dummy_atom;
    } else {
        frag_other_atom = get_atom_bound_to_ap_parent_atom(
            fragment, frag_ap_parent_atom, is_ring);
    }

    RDKit::Conformer frag_conf =
        translate_fragment_ap_to(fragment, core_atom_coords);
    auto frag_other_atom_coords =
        frag_conf.getAtomPos(frag_other_atom->getIdx());
    auto rotate_by = get_angle_radians(frag_other_atom_coords, core_atom_coords,
                                       core_best_placing);
    if (is_ring) {
        // rotate_by will put the attachment point in the best position, instead
        // of the ring itself, so we flip the angle 180 degrees
        rotate_by += M_PI;
    }
    rotate_conformer_radians(rotate_by, core_atom_coords, frag_conf);
    return frag_conf;
}

std::pair<RDKit::Conformer, const RDKit::Atom*>
align_fragment_with_bond(const RDKit::ROMol& fragment,
                         const RDKit::Bond* const core_bond)
{
    const auto& core = core_bond->getOwningMol();
    auto [frag_ap_parent_atom, frag_ap_dummy_atom, is_ring] =
        get_frag_info(fragment);
    const RDKit::Bond* frag_bond = get_frag_bond_to_align_to(
        fragment, frag_ap_parent_atom, core, core_bond);
    const RDKit::Atom* core_atom; // the core atom to overlay frag_atom on
    if (is_ring) {
        // Arbitrarily pick one of core_bond's atoms to overlay.  We'll test the
        // opposite configuration (i.e. rotating 180 degrees) below.
        core_atom = core_bond->getBeginAtom();
    } else {
        // use the atom with the highest degree as the point to overlay so that
        // the fragment itself will always be over any terminal atoms if the
        // user hovers over a terminal bond
        core_atom = std::max(core_bond->getBeginAtom(), core_bond->getEndAtom(),
                             [&core](const RDKit::Atom* const atom1,
                                     const RDKit::Atom* const atom2) {
                                 return core.getAtomDegree(atom1) <
                                        core.getAtomDegree(atom2);
                             });
    }

    auto new_conf = overlay_fragment_on_core(
        fragment, frag_ap_parent_atom, frag_bond, core, core_atom, core_bond);
    if (is_ring &&
        should_flip_ring_fragment(fragment, new_conf, frag_ap_dummy_atom, core,
                                  core_bond)) {
        // flip new_conf by 180 degrees so that the ring points away from
        // neighboring core bonds
        core_atom = core_bond->getOtherAtom(core_atom);
        new_conf =
            overlay_fragment_on_core(fragment, frag_ap_parent_atom, frag_bond,
                                     core, core_atom, core_bond);
    }
    return {new_conf, core_atom};
}

std::unordered_set<RDKit::Bond*>
determine_fragment_bonds_to_mutate_to_single_bonds(
    RDKit::ROMol& fragment, const RDKit::Conformer& frag_conf,
    const RDKit::ROMol& core, const RDKit::Atom* const core_start_atom)
{
    auto [frag_ap_parent_atom, frag_ap_dummy_atom, is_ring] =
        get_frag_info(fragment);
    auto frag_atom_to_core_atom = get_fragment_to_core_atom_map(
        fragment, frag_conf, frag_ap_parent_atom, core, core_start_atom);
    return determine_fragment_bonds_to_mutate_to_single_bonds(
        fragment, core, frag_atom_to_core_atom);
}

/**
 * @overload An overload that takes a mapping between core and fragment atoms
 * instead of using the conformers.
 */
static std::unordered_set<RDKit::Bond*>
determine_fragment_bonds_to_mutate_to_single_bonds(
    RDKit::ROMol& fragment, const RDKit::ROMol& core,
    const AtomToAtomMap& frag_atom_to_core_atom)
{
    std::unordered_set<RDKit::Bond*> frag_bonds_to_mutate;
    auto bond_order_gt_one = [](const auto* bond) {
        return bond->getBondTypeAsDouble() > 1;
    };
    const auto* table = RDKit::PeriodicTable::getTable();
    for (auto [cur_frag_atom, cur_core_atom] : frag_atom_to_core_atom) {
        if (fragment.getAtomDegree(cur_frag_atom) <= 1) {
            continue;
        }
        auto cur_final_atom =
            should_replace_core_atom(cur_frag_atom, cur_core_atom)
                ? cur_frag_atom
                : cur_core_atom;
        auto valences = table->getValenceList(cur_final_atom->getAtomicNum());
        auto max_valence = *std::max_element(valences.begin(), valences.end());

        float cur_valence = 0.0;
        std::unordered_set<const RDKit::Bond*>
            core_bonds_that_overlap_frag_bonds;

        // count the total valence contribution of the fragment bonds (and any
        // core bonds that overlap fragment bonds) at this atom
        for (auto* frag_bond : fragment.atomBonds(cur_frag_atom)) {
            auto [frag_neighbor, core_neighbor, core_bond] =
                get_frag_neighbor_and_corresponding_core_neighbor_and_bond(
                    cur_frag_atom, frag_bond, cur_core_atom,
                    frag_atom_to_core_atom);
            if (core_bond == nullptr ||
                should_replace_core_bond(frag_bond, core_bond)) {
                cur_valence += frag_bond->getBondTypeAsDouble();
            } else {
                cur_valence += core_bond->getBondTypeAsDouble();
            }
            if (core_bond != nullptr) {
                core_bonds_that_overlap_frag_bonds.insert(core_bond);
                if (core_bond->getBondType() == RDKit::Bond::BondType::DOUBLE &&
                    frag_bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
                    // Mutate this fragment bond to avoid showing a blue
                    // fragment hint double bond on top of a core double bond.
                    // If the bonds are asymmetrical, then they'll wind up
                    // looking like a triple bond when they're overlayed, since
                    // the fragment bond's second line will be on the opposite
                    // side from the core bond's second line. Note that this
                    // mutation will *not* affect the actual structure that's
                    // inserted, since the core's double bond will "win" over
                    // the fragment's (post-mutation) single bond.
                    frag_bonds_to_mutate.insert(frag_bond);
                }
            }
        }

        // count the total valence contribution of the remaining core bonds at
        // this atom
        for (auto* core_bond : core.atomBonds(cur_core_atom)) {
            if (!core_bonds_that_overlap_frag_bonds.count(core_bond)) {
                cur_valence += core_bond->getBondTypeAsDouble();
            }
        }

        if (cur_valence > max_valence) {
            // this atom is going to have too high of a valence if we don't do
            // any mutations, so mutate all non-single fragment bonds to single
            // bonds
            auto frag_bonds_it = fragment.atomBonds(cur_frag_atom);
            std::copy_if(frag_bonds_it.begin(), frag_bonds_it.end(),
                         std::inserter(frag_bonds_to_mutate,
                                       frag_bonds_to_mutate.begin()),
                         bond_order_gt_one);
        }
    }
    return frag_bonds_to_mutate;
}

/**
 * @brief Find the core bond, if any, that corresponds to the given fragment
 * bond
 * @param frag_atom The fragment atom that corresponds to core_atom
 * @param frag_bond A fragment bond involving frag_atom
 * @param core_atom The core atom that corresponds to frag_atom
 * @param frag_atom_to_core_atom  A map of fragment atoms to their
 * corresponding core atom
 * @return A tuple of
 *   - The other fragment atom (i.e. not frag_atom) from frag_bond
 *   - The core atom that corresponds to the other fragment atom from frag_bond.
 *     Will be nullptr if there is no corresponding atom.
 *   - The core bond that corresponds to frag_bond. Will be nullptr if there is
 *     no corresponding bond.
 */
static std::tuple<const RDKit::Atom*, const RDKit::Atom*, const RDKit::Bond*>
get_frag_neighbor_and_corresponding_core_neighbor_and_bond(
    const RDKit::Atom* const frag_atom, const RDKit::Bond* const frag_bond,
    const RDKit::Atom* const core_atom,
    const AtomToAtomMap& frag_atom_to_core_atom)
{
    auto* frag_neighbor = frag_bond->getOtherAtom(frag_atom);
    const RDKit::Atom* core_neighbor = nullptr;
    const RDKit::Bond* core_bond = nullptr;
    if (frag_atom_to_core_atom.count(frag_neighbor)) {
        core_neighbor = frag_atom_to_core_atom.at(frag_neighbor);
        auto& core = core_atom->getOwningMol();
        core_bond = core.getBondBetweenAtoms(core_atom->getIdx(),
                                             core_neighbor->getIdx());
    }
    return {frag_neighbor, core_neighbor, core_bond};
}

/**
 * Determine whether we should replace a core bond with the corresponding
 * fragment bond
 * @param frag_bond The fragment bond that corresponds to core_bond
 * @param core_bond The core bond that corresponds to frag_bond
 */
static bool should_replace_core_bond(const RDKit::Bond* const frag_bond,
                                     const RDKit::Bond* const core_bond)
{
    return frag_bond->getBondType() != RDKit::Bond::BondType::SINGLE;
}

RDKit::Atom* prepare_fragment_for_insertion(RDKit::RWMol& fragment)
{
    auto [frag_ap_parent_atom, frag_ap_dummy_atom, is_ring] =
        get_frag_info(fragment);
    // remove the attachment point atom, which will automatically remove the
    // attachment point bond
    fragment.removeAtom(frag_ap_dummy_atom->getIdx());
    update_molecule_on_change(fragment);
    // strip constness from frag_ap_parent_atom, since get_frag_info returns
    // const values
    return fragment.getAtomWithIdx(frag_ap_parent_atom->getIdx());
}

/**
 * When adding a fragment, figure out which fragment atoms correspond to core
 * atoms.  In order to "correspond," the two atoms must have the same
 * coordinates.  Additionally, all corresponding atoms must be connected.  In
 * other words, frag_start_atom always corresponds to core_start_atom.
 * Neighbors of frag_start_atom may correspond to neighbors of core_start_atom,
 * assuming they have the same coordinates.  If a pair of those atoms are
 * corresponding, then any neighbors of *those* atoms may correspond if they
 * have the same coordinates, etc.
 * @param fragment The fragment being added
 * @param frag_conf The conformer of the fragment
 * @param frag_start_atom The "starting" atom for the fragment, i.e. the
 * attachment point parent atom.  This atom must have the same coordinates as
 * core_start_atom.
 * @param core The core molecule
 * @param core_start_atom The "starting" atom for the core, i.e. the atom
 * that was clicked on.  This atom must have the same coordinates as
 * frag_start_atom.
 * @return A map of fragment atoms to their corresponding core atom
 */
static AtomToAtomMap get_fragment_to_core_atom_map(
    const RDKit::ROMol& fragment, const RDKit::Conformer& frag_conf,
    const RDKit::Atom* frag_start_atom, const RDKit::ROMol& core,
    const RDKit::Atom* core_start_atom)
{
    AtomToAtomMap frag_atom_to_core_atom{{frag_start_atom, core_start_atom}};
    std::queue<std::pair<const RDKit::Atom*, const RDKit::Atom*>> to_examine;
    to_examine.emplace(frag_start_atom, core_start_atom);

    const auto& core_conf = core.getConformer();
    while (!to_examine.empty()) {
        auto [frag_atom, core_atom] = to_examine.front();
        to_examine.pop();
        for (auto* frag_neighbor : fragment.atomNeighbors(frag_atom)) {
            if (frag_atom_to_core_atom.count(frag_neighbor)) {
                // we've already examined this atom and found a match, so
                // there's nothing more to do
                continue;
            }
            const auto& frag_neighbor_coords =
                frag_conf.getAtomPos(frag_neighbor->getIdx());
            for (auto* core_neighbor : core.atomNeighbors(core_atom)) {
                const auto& core_neighbor_coords =
                    core_conf.getAtomPos(core_neighbor->getIdx());
                auto dist =
                    (frag_neighbor_coords - core_neighbor_coords).length();
                if (dist <= MAX_DIST_FOR_FRAG_OVERLAP) {
                    frag_atom_to_core_atom[frag_neighbor] = core_neighbor;
                    to_examine.emplace(frag_neighbor, core_neighbor);
                }
            }
        }
    }
    return frag_atom_to_core_atom;
}

/**
 * @overload An overload of that uses the default conformer for the fragment
 */
static AtomToAtomMap get_fragment_to_core_atom_map(
    const RDKit::ROMol& fragment, const RDKit::Atom* frag_start_atom,
    const RDKit::ROMol& core, const RDKit::Atom* core_start_atom)
{
    const auto& frag_conf = fragment.getConformer();
    return get_fragment_to_core_atom_map(fragment, frag_conf, frag_start_atom,
                                         core, core_start_atom);
}

/**
 * Determine whether we should replace a core atom with a fragment atom
 * @param frag_atom The fragment atom that corresponds to core_atom
 * @param core_atom The core atom that corresponds to frag_atom
 */
static bool should_replace_core_atom(const RDKit::Atom* const frag_atom,
                                     const RDKit::Atom* const core_atom)
{
    if (is_heteroatom_or_query(frag_atom)) {
        return true;
    }
    return !is_heteroatom_or_query(core_atom);
}

/**
 * When adding a fragment, figure out which core atoms should be mutated to
 * match the fragment.
 * @param frag_atom_to_core_atom A map of fragment atoms to their
 * corresponding core atom
 * @return A list of (core atom index, function that returns an atom to mutate
 * to)
 */
static std::vector<std::pair<unsigned int, AtomFunc>>
determine_core_atom_mutations(const AtomToAtomMap& frag_atom_to_core_atom)
{
    std::vector<std::pair<unsigned int, AtomFunc>> mutations_to_core_atoms;
    for (auto [frag_atom, core_atom] : frag_atom_to_core_atom) {
        if (should_replace_core_atom(frag_atom, core_atom)) {
            auto& frag_atom_ref = *frag_atom;
            auto copy_frag_atom = [frag_atom_ref]() {
                return std::make_shared<RDKit::Atom>(frag_atom_ref);
            };
            mutations_to_core_atoms.emplace_back(core_atom->getIdx(),
                                                 copy_frag_atom);
        }
    }
    return mutations_to_core_atoms;
}

/**
 * When adding a fragment, figure out which bonds involving at least one
 * core atom should be added or mutated to match the fragment.
 * @param core The core molecule
 * @param fragment The fragment structure being added
 * @param frag_atom_to_core_atom A map of fragment atoms to their
 * corresponding core atom
 * @return A tuple of:
 *   - Bonds to make between the core and the fragment, formatted as a map
 *     of {core atom: {fragment atom: info about bond to create}}
 *   - A list of core bonds to mutate, where each bond is given as (core
 *     bond index, function that returns a bond to mutate to)
 *   - A list of bonds to make between two core atoms, where each bond is
 *     given as a tuple of
 *     - atom index for the starting atom
 *     - atom index for the ending atom
 *     - function that returns a bond instance for the new bond
 */
static std::tuple<AtomPtrToFragBondMap,
                  std::vector<std::pair<unsigned int, BondFunc>>,
                  std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>>
determine_core_bond_additions_and_mutations(
    const RDKit::ROMol& core, const RDKit::ROMol& fragment,
    AtomToAtomMap frag_atom_to_core_atom)
{
    AtomPtrToFragBondMap core_to_frag_bonds;
    std::vector<std::pair<unsigned int, BondFunc>> mutations_to_core_bonds;
    std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>
        additions_to_core_bonds;
    for (auto [frag_atom, core_atom] : frag_atom_to_core_atom) {
        for (auto* frag_bond : fragment.atomBonds(frag_atom)) {
            auto [frag_neighbor, core_neighbor, core_bond] =
                get_frag_neighbor_and_corresponding_core_neighbor_and_bond(
                    frag_atom, frag_bond, core_atom, frag_atom_to_core_atom);
            auto& frag_bond_ref = *frag_bond;
            auto copy_bond = [frag_bond_ref]() {
                return std::make_shared<RDKit::Bond>(frag_bond_ref);
            };
            if (core_neighbor != nullptr) {
                // both of the fragment atoms in this bond are on top of a core
                // atom
                if (core_bond == nullptr) {
                    // there's no equivalent bond in the core, so one needs to
                    // be added
                    int start_atom_idx = core_atom->getIdx();
                    int end_atom_idx = core_neighbor->getIdx();
                    if (frag_bond->getBeginAtom() != frag_atom) {
                        std::swap(start_atom_idx, end_atom_idx);
                    }
                    additions_to_core_bonds.emplace_back(
                        start_atom_idx, end_atom_idx, copy_bond);
                } else if (should_replace_core_bond(frag_bond, core_bond)) {
                    // there's already a core bond between these atoms but the
                    // fragment bond is "interesting" (i.e. not a single bond),
                    // so we mutate the core bond to match the fragment bond
                    mutations_to_core_bonds.emplace_back(core_bond->getIdx(),
                                                         copy_bond);
                }
            } else {
                // this bond will be between a core atom and a fragment atom
                core_to_frag_bonds.try_emplace(core_atom);
                bool frag_first = frag_bond->getBeginAtom() == frag_atom;
                FragmentBondInfo bond_info =
                    std::make_tuple(copy_bond, frag_first);
                core_to_frag_bonds[core_atom].emplace(frag_neighbor, bond_info);
            }
        }
    }
    return {core_to_frag_bonds, mutations_to_core_bonds,
            additions_to_core_bonds};
}

/**
 * Convert a map of atoms-to-bond-info from using atom pointers to using atom
 * indices.  For the fragment, the atom indices used are the indices that the
 * fragment atoms *will have* once they are added to the core molecule.
 * @param core The core molecule
 * @param core_to_frag_bonds_by_ptr A map of {core atom: {fragment atom:
 * info about bond to create}}, where both atoms are represented by pointers
 * @return A map that is identical to core_to_frag_bonds_by_ptr, but using
 * atom indices in place of pointers.
 */
static AtomIdxToFragBondMap convert_bond_map_from_ptrs_to_idxs(
    const RDKit::ROMol& core,
    const AtomPtrToFragBondMap& core_to_frag_bonds_by_ptr)
{
    AtomIdxToFragBondMap core_to_frag_bonds_by_idx;
    unsigned int num_core_atoms = core.getNumAtoms();
    for (auto& [core_atom, frag_atom_to_bond_info] :
         core_to_frag_bonds_by_ptr) {
        auto core_atom_idx = core_atom->getIdx();
        core_to_frag_bonds_by_idx.try_emplace(core_atom_idx);
        auto& frag_atom_idx_to_bond_info =
            core_to_frag_bonds_by_idx[core_atom_idx];
        for (auto& [frag_atom, bond_info] : frag_atom_to_bond_info) {
            auto frag_atom_idx = num_core_atoms + frag_atom->getIdx();
            frag_atom_idx_to_bond_info.emplace(frag_atom_idx, bond_info);
        }
    }
    return core_to_frag_bonds_by_idx;
}

std::tuple<std::vector<std::pair<unsigned int, AtomFunc>>,
           std::vector<std::pair<unsigned int, BondFunc>>,
           std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>,
           AtomIdxToFragBondMap>
get_fragment_addition_info(RDKit::RWMol& fragment,
                           const RDKit::Atom* frag_start_atom,
                           const RDKit::ROMol& core,
                           const RDKit::Atom* core_start_atom)
{
    // figure out which fragment atoms correspond to which atoms in the core
    // (i.e. the existing structure)
    auto frag_atom_to_core_atom = get_fragment_to_core_atom_map(
        fragment, frag_start_atom, core, core_start_atom);

    // convert fragment bonds to single bonds to try to avoid valence errors on
    // the core atoms where the fragment is added
    auto frag_bonds_to_mutate =
        determine_fragment_bonds_to_mutate_to_single_bonds(
            fragment, core, frag_atom_to_core_atom);
    for (auto* cur_bond : frag_bonds_to_mutate) {
        cur_bond->setBondType(RDKit::Bond::BondType::SINGLE);
    }

    // figure out which core atoms should be mutated so they match the overlayed
    // fragment atom
    auto mutations_to_core_atoms =
        determine_core_atom_mutations(frag_atom_to_core_atom);

    // figure out where bonds should be added and which core bonds should be
    // mutated
    auto [core_to_frag_bonds_by_ptr, mutations_to_core_bonds,
          additions_to_core_bonds] =
        determine_core_bond_additions_and_mutations(core, fragment,
                                                    frag_atom_to_core_atom);

    // delete all fragment atoms that overlay on existing core atoms
    for (auto [frag_atom, core_atom] : frag_atom_to_core_atom) {
        fragment.removeAtom(frag_atom->getIdx());
    }

    // core_to_frag_bonds_by_ptr uses atom pointers, so it's only valid for this
    // exact instance of core and fragment.  Convert the pointers to atom
    // indices so that the dictionary will still be valid even when using a copy
    // of core or frag (which happens in MolModel::addFragmentFromCommand).
    auto core_to_frag_bonds_by_idx =
        convert_bond_map_from_ptrs_to_idxs(core, core_to_frag_bonds_by_ptr);

    return {mutations_to_core_atoms, mutations_to_core_bonds,
            additions_to_core_bonds, core_to_frag_bonds_by_idx};
}

} // namespace sketcher
} // namespace schrodinger
