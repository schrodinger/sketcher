/// This file contains utility functions used in more than one translation unit
/// related to monomer conversion etc., but which are not meant as part of the
/// public (Python) API.

#pragma once

#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/definitions.h"

#include <rdkit/GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

inline constexpr unsigned int NO_ATTACHMENT =
    std::numeric_limits<unsigned int>::max();

/// Return the atom map number for each atom in the molecule, with a default
/// of NO_ATTACHMENT for atoms without a map number.
[[nodiscard]] std::vector<unsigned int> make_attch_map(const RDKit::ROMol& mol);

void neutralizeAtoms(RDKit::ROMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
