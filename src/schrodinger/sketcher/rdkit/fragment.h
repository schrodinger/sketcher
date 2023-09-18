/**
 * Functions for working with fragments used by DrawFragmentSceneTool (e.g. the
 * standard draw ring tools).  Note that all fragments must have exactly one
 * attachment point, as determined by the is_attachment_point method in
 * rgroup.h.
 */

#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

typedef std::function<std::shared_ptr<RDKit::Atom>()> AtomFunc;
typedef std::function<std::shared_ptr<RDKit::Bond>()> BondFunc;

typedef std::unordered_map<const RDKit::Atom*, const RDKit::Atom*>
    AtomToAtomMap;

// a tuple of:
//   - function that returns a Bond instance to use for the bond
//   - whether the fragment atom is the starting atom of the bond (true) or the
//     ending atom (false)
typedef std::tuple<BondFunc, bool> FragmentBondInfo;

typedef std::unordered_map<
    const RDKit::Atom*,
    std::unordered_map<const RDKit::Atom*, FragmentBondInfo>>
    AtomPtrToFragBondMap;
typedef std::unordered_map<unsigned int,
                           std::unordered_map<unsigned int, FragmentBondInfo>>
    AtomIdxToFragBondMap;

/**
 * Translate the fragment so that the attachment point parent atom (i.e. the
 * heavy atom bound to the attachment point dummy atom) is at the specified
 * point
 * @return A newly generated conformer.  Note that this conformer has *not* been
 * added to the molecule.
 */
SKETCHER_API RDKit::Conformer translate_fragment(const RDKit::ROMol& fragment,
                                                 const RDGeom::Point3D& point);

/**
 * Align the fragment so that it grows out of the specified atom
 * @param fragment The fragment to align
 * @param core_atom The atom to use.  Note that this atom must not be an
 * attachment point and must not be part of fragment.
 * @return A new fragment conformer.  Note that this conformer is not added to
 * fragment.
 */
SKETCHER_API RDKit::Conformer
align_fragment_with_atom(const RDKit::ROMol& fragment,
                         const RDKit::Atom* const core_atom);

/**
 * Align the fragment so that it grows out of the specified bond
 * @param fragment The fragment to align
 * @param core_bond The bond to use.  Note that this bond must not be an
 * attachment point bond and must not be part of fragment.
 * @return A new fragment conformer.  Note that this conformer is not added to
 * fragment.
 */
SKETCHER_API std::pair<RDKit::Conformer, const RDKit::Atom*>
align_fragment_with_bond(const RDKit::ROMol& fragment,
                         const RDKit::Bond* const core_bond);

/**
 * Prepare the given fragment for insertion into the core molecule by removing
 * the attachment point and attachment point bond.
 * @param[in,out] fragment The fragment to prepare
 * @return The fragment's attachment point parent atom (i.e. the atom that was
 * bound to the now-removed attachment point dummy atom).  This atom should
 * correspond to the clicked core atom.
 */
SKETCHER_API RDKit::Atom*
prepare_fragment_for_insertion(RDKit::RWMol& fragment);

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
 * @param frag_start_atom The "starting" atom for the fragment, i.e. the
 * attachment point parent atom.  This atom must have the same coordinates as
 * core_start_atom.
 * @param core The core molecule
 * @param core_start_atom The "starting" atom for the core, i.e. the atom
 * that was clicked on.  This atom must have the same coordinates as
 * frag_start_atom.
 * @return A map of fragment atoms to their corresponding core atom
 */
SKETCHER_API AtomToAtomMap get_fragment_to_core_atom_map(
    const RDKit::ROMol& fragment, const RDKit::Atom* frag_start_atom,
    const RDKit::ROMol& core, const RDKit::Atom* core_start_atom);

/**
 * Determine whether we should replace a core atom with a fragment atom
 * @param frag_atom The fragment atom that corresponds to core_atom
 * @param core_atom The core atom that corresponds to frag_atom
 */
SKETCHER_API bool should_replace_core_atom(const RDKit::Atom* const frag_atom,
                                           const RDKit::Atom* const core_atom);

/**
 * When adding a fragment, figure out which core atoms should be mutated to
 * match the fragment.
 * @param frag_atom_to_core_atom A map of fragment atoms to their
 * corresponding core atom
 * @return A list of (core atom index, function that returns an atom to mutate
 * to)
 */
std::vector<std::pair<unsigned int, AtomFunc>>
determine_core_atom_mutations(const AtomToAtomMap& frag_atom_to_core_atom);

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
std::tuple<AtomPtrToFragBondMap, std::vector<std::pair<unsigned int, BondFunc>>,
           std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>>
determine_core_bond_additions_and_mutations(
    const RDKit::ROMol& core, const RDKit::ROMol& fragment,
    AtomToAtomMap frag_atom_to_core_atom);

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
AtomIdxToFragBondMap convert_bond_map_from_ptrs_to_idxs(
    const RDKit::ROMol& core,
    const AtomPtrToFragBondMap& core_to_frag_bonds_by_ptr);

/**
 * Figure out all changes that need to be made in order to connect a fragment to
 * an existing core molecule.
 * @param[in,out] fragment The fragment being added.  This function will delete
 * any atoms from this molecule that should be replaced by the corresponding
 * core atom.  See the get_fragment_to_core_atom_map docstring for an
 * explanation of corresponding atoms.
 * @param[in] frag_start_atom The "starting" atom for the fragment, i.e. the
 * attachment point parent atom.  This atom must have the same coordinates as
 * core_start_atom.
 * @param[in] core The core molecule
 * @param[in] core_start_atom The "starting" atom for the core, i.e. the atom
 * that was clicked on.  This atom must have the same coordinates as
 * frag_start_atom.
 */
std::tuple<std::vector<std::pair<unsigned int, AtomFunc>>,
           std::vector<std::pair<unsigned int, BondFunc>>,
           std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>,
           AtomIdxToFragBondMap>
get_fragment_addition_info(RDKit::RWMol& fragment,
                           const RDKit::Atom* frag_start_atom,
                           const RDKit::ROMol& core,
                           const RDKit::Atom* core_start_atom);

} // namespace sketcher
} // namespace schrodinger
