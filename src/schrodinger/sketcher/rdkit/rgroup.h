#pragma once

#include <memory>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
class Bond;
class Conformer;
class ROMol;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a new attachment point atom using the specified attachment point
 * number
 */
std::shared_ptr<RDKit::Atom>
make_new_attachment_point(const unsigned int ap_num);

/**
 * @return whether the specified atom is an attachment point.  Note that we
 * require attachment points to have an assigned number, so the atom label must
 * be "_AP" followed by a number, not just "_AP".
 */
SKETCHER_API bool is_attachment_point(const RDKit::Atom* const atom);

/**
 * @return The attachment point number of the specified atom.  If the atom is
 * not an attachment point, 0 will be returned.
 */
SKETCHER_API unsigned int
get_attachment_point_number(const RDKit::Atom* const atom);

/**
 * @return whether the specified atom is an R-group.  Note that we require
 * R-groups to have an assigned number, so the atom label must be "_R" followed
 * by a number.
 */
SKETCHER_API bool is_r_group(const RDKit::Atom* const atom);

/**
 * @return A sorted list of all used R-group numbers in the molecule
 */
SKETCHER_API std::vector<unsigned int>
get_all_r_group_numbers(const RDKit::ROMol* const mol);

/**
 * Get multiple available R-group numbers
 * @param how_many How many R-group numbers to return
 * @return A sorted list (of length `how_many`) containing the smallest
 * unused R-group numbers
 */
SKETCHER_API std::vector<unsigned int>
get_next_r_group_numbers(const RDKit::ROMol* const mol, const size_t how_many);

/**
 * @return The next available (i.e. smallest unused) attachment point number
 */
SKETCHER_API unsigned int
get_next_attachment_point_number(const RDKit::ROMol* const mol);

/**
 * Number all access points consecutively
 */
SKETCHER_API void renumber_attachment_points(RDKit::RWMol* const mol);

/**
 * @return the bond bound to a given attachment point dummy atom.  If a
 * non-attachment point atom is passed in, nullptr will be returned.
 */
SKETCHER_API const RDKit::Bond* const
get_attachment_point_bond(const RDKit::Atom* const atom);

/**
 * @return the attachment point dummy atom bound to the given bond.  If a
 * non-attachment point bond is passed in, nullptr will be returned.
 */
SKETCHER_API const RDKit::Atom* const
get_attachment_point_atom(const RDKit::Bond* const bond);

/**
 * @return whether the given bond is bound to an attachment point dummy atom.
 */
SKETCHER_API bool is_attachment_point_bond(const RDKit::Bond* const bond);

/**
 * @return The number of attachment point dummy atoms bound to the specified
 * atom
 */
SKETCHER_API unsigned int
number_of_bound_attachment_points(const RDKit::Atom* const atom);

/**
 * Shorten all bonds bound to attachment point dummy atoms
 * @param mol The conformer to modify.
 */
SKETCHER_API void shorten_attachment_point_bonds(RDKit::Conformer& conf);

} // namespace sketcher
} // namespace schrodinger
