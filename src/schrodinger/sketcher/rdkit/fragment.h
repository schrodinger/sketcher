#include <GraphMol/ROMol.h>

namespace schrodinger
{
namespace sketcher
{

// TODO: report whether the AP is part of a ring?
/**
 * Get the atom bound to the attachment point
 * @param fragment A molecule with exactly one attachment point
 */
RDKit::Atom* const
get_attachment_point_heavy_atom(const RDKit::ROMol& fragment);

/**
 * Translate the fragment so that the attachment point heavy atom is at the
 * specified point
 * @return A newly generated conformer.  Note that this conformer has *not* been
 * added to the molecule.
 */
RDKit::Conformer translate_fragment(const RDKit::ROMol& fragment,
                                    const RDGeom::Point3D& point);

// TODO
// TODO: report overlapping atoms?
RDKit::Conformer align_fragment_with_atom(const RDKit::ROMol& fragment,
                                          const RDKit::Atom* const atom);

// TODO
// TODO: report overlapping atoms?
RDKit::Conformer align_fragment_with_bond(const RDKit::ROMol& fragment,
                                          const RDKit::Bond* const bond);

} // namespace sketcher
} // namespace schrodinger
