/**
 * Functions for working with fragments used by DrawFragmentSceneTool (e.g. the
 * standard draw ring tools).  Note that all fragments must have exactly one
 * attachment point, as determined by the is_attachment_point method in
 * rgroup.h.
 */

#include <rdkit/GraphMol/ROMol.h>

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
SKETCHER_API RDKit::Conformer
translate_fragment_ap_to(const RDKit::ROMol& fragment,
                         const RDGeom::Point3D& point);

/**
 * Translate the fragment so that the fragment center (calculated over all atoms
 * *other than* the attachment point dummy atom) is at the specified point
 * @return A newly generated conformer.  Note that this conformer has *not* been
 * added to the molecule.
 */
SKETCHER_API RDKit::Conformer
translate_fragment_center_to(const RDKit::ROMol& fragment,
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
 * Determine which fragment bonds should be mutated to single bonds to avoid
 * valence errors where the fragment attaches to the core.  For each fragment
 * atom that overlays on a core atom, we convert all of the fragment atom's
 * bonds to single bonds if both of the following are true:
 *   - The core atom that was clicked on has an aromatic, double, triple, or
 *     quadruple bond
 *   - The fragment atom has more than one bond, at least one of which is an
 *     aromatic, double, triple, or quadruple bond
 * @param fragment The fragment structure to attach
 * @param frag_conf The conformer that positions the fragment as it will be
 * attached
 * @param core The core structure
 * @param core_start_atom A core atom where the fragment will be attached.  This
 * is typically the core atom that was clicked on, or one of the atoms in the
 * bond that was clicked on.
 * @return A list of fragments bonds to mutate
 */
SKETCHER_API std::unordered_set<RDKit::Bond*>
determine_fragment_bonds_to_mutate_to_single_bonds(
    RDKit::ROMol& fragment, const RDKit::Conformer& frag_conf,
    const RDKit::ROMol& core, const RDKit::Atom* const core_start_atom);

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
