#pragma once

#include <optional>
#include <string>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

// Forward declarations:
namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class EnhancedStereo;

/**
 * @param strip_abs if true, the absolute stereo label is stripped from the
 * "abs" label
 * @param show_unspecified if true, a (?) label is shown on unspecified
 * stereo centers
 * @return the chiral label for the given atom.
 */
SKETCHER_API std::string get_atom_chirality_label(const RDKit::Atom& atom,
                                                  bool strip_abs = false,
                                                  bool show_unspecified = true);

/**
 * @return the chiral label for the given bond.
 */
SKETCHER_API std::string get_bond_stereo_label(const RDKit::Bond& bond);

/**
 * Get the enhanced stereo value for the specified atom
 * @return The enhanced stereo type and group write ID for the given atom. (Note
 * that the group ID is only valid for AND or OR stereochemistry, as absolute
 * stereochemistry does not have a group id.) If the atom is not part of any
 * enhanced stereo group, then std::nullopt will be returned.
 */
SKETCHER_API std::optional<EnhancedStereo>
get_enhanced_stereo_for_atom(const RDKit::Atom* atom);

/**
 * Modify the enhanced stereo group for the specified atom
 * @param atom The atom to modify
 * @param enh_stereo A pair of the enhanced stereo type and group ID. For AND
 * and OR stereochemistry, this group ID will be set for both the read and write
 * ID. For absolute stereochemistry, this group ID will be ignored.
 */
SKETCHER_API void
set_enhanced_stereo_for_atom(RDKit::Atom* atom,
                             const EnhancedStereo& enh_stereo);

} // namespace sketcher
} // namespace schrodinger
