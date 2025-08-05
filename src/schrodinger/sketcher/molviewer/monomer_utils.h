/**
 * Utility functions for flagging monomeric atoms in MolModel. Note that these
 * functions *are not generally relevant* to RDKit molecules. (For example,
 * is_atom_monomeric() only looks for the flag set by set_atom_monomeric(), and
 * that flag won't be set by anything outside of MolModel.)
 */

#pragma once

#include <string>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Flag the given atom as monomeric
 */
SKETCHER_API void set_atom_monomeric(RDKit::Atom* const atom);

/**
 * Remove the monomeric Flag from the atom if one is present. This should be
 * called before exporting the mol to avoid leaking MolModel internal flags.
 */
SKETCHER_API void clear_monomeric_property(RDKit::Atom* const atom);

/**
 * @return whether the atom has been flagged as monomeric (via
 * set_atom_monomeric)
 */
SKETCHER_API bool is_atom_monomeric(const RDKit::Atom* const atom);

/**
 * @return whether any atoms in the mol have been flagged as monomeric (via
 * set_atom_monomeric)
 */
SKETCHER_API bool contains_monomeric_atom(const RDKit::ROMol& mol);

} // namespace sketcher
} // namespace schrodinger
