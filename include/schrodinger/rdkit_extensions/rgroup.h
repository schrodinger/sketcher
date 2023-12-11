/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: rgroup querying and creations apis
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>
#include <optional>

#include "schrodinger/rdkit_extensions/definitions.h" // RDKIT_EXTENSIONS_API

// Forward declarations:
namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @param atom rdkit atom
 * @return whether the atom is an attachment point
 */
RDKIT_EXTENSIONS_API bool is_attachment_point_dummy(const RDKit::Atom& atom);

/**
 * @return The R-group number of the specified atom.  If the atom is not an
 * R-group, an empty `r_group_num_t` object will be returned.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<unsigned int>
get_r_group_number(const RDKit::Atom* const atom);

/**
 * @return a new R-group atom using the specified R-group number
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::shared_ptr<RDKit::Atom>
make_new_r_group(const unsigned int r_group_num);

} // namespace rdkit_extensions
} // namespace schrodinger
