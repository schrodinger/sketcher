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

typedef std::unordered_map<const RDKit::Atom*, const RDKit::Atom*>
    AtomToAtomMap;

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

} // namespace sketcher
} // namespace schrodinger
