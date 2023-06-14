#pragma once

#include <memory>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class ROMol;
class RWMol;
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a new R-group atom using the specified R-group number
 */
std::shared_ptr<RDKit::Atom> make_new_r_group(const unsigned int r_group_num);

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
 * @return The R-group number of the specified atom.  If the atom is not an
 * R-group, 0 will be returned.
 */
SKETCHER_API unsigned int get_r_group_number(const RDKit::Atom* const atom);

/**
 * @return A sorted list of all used R-group numbers in the molecule
 */
SKETCHER_API std::vector<unsigned int>
get_all_r_group_numbers(const RDKit::ROMol* const mol);

/**
 * @return The next available (i.e. smallest unused) R-group number
 */
SKETCHER_API unsigned int
get_next_r_group_number(const RDKit::ROMol* const mol);

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

} // namespace sketcher
} // namespace schrodinger
