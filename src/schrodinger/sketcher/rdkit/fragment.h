/**
 * Functions for working with fragments used by DrawFragmentSceneTool (e.g. the
 * standard draw ring tools).  Note that all fragments must have exactly one
 * attachment point, as determined by the is_attachment_point method in
 * rgroup.h.
 */

#include <GraphMol/ROMol.h>

namespace schrodinger
{
namespace sketcher
{

/**
 * Translate the fragment so that the attachment point parent atom (i.e. the
 * heavy atom bound to the attachment point dummy atom) is at the specified
 * point
 * @return A newly generated conformer.  Note that this conformer has *not* been
 * added to the molecule.
 */
RDKit::Conformer translate_fragment(const RDKit::ROMol& fragment,
                                    const RDGeom::Point3D& point);

/**
 * Align the fragment so that it grows out of the specified atom
 * @param fragment The fragment to align
 * @param core_atom The atom to use.  Note that this atom must not be an
 * attachment point and must not be part of fragment.
 * @return A new fragment conformer.  Note that this conformer is not added to
 * fragment.
 */
RDKit::Conformer align_fragment_with_atom(const RDKit::ROMol& fragment,
                                          const RDKit::Atom* const core_atom);

/**
 * Align the fragment so that it grows out of the specified bond
 * @param fragment The fragment to align
 * @param core_bond The bond to use.  Note that this bond must not be an
 * attachment point bond and must not be part of fragment.
 * @return A new fragment conformer.  Note that this conformer is not added to
 * fragment.
 */

RDKit::Conformer align_fragment_with_bond(const RDKit::ROMol& fragment,
                                          const RDKit::Bond* const core_bond);

} // namespace sketcher
} // namespace schrodinger
