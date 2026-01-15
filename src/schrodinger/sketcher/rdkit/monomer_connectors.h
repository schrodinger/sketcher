#pragma once

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * @return whether the given bond represents two connections between the same
 * monomers, such as neighboring cysteines additionally joined by a disulfide
 * bond. RDKit does not allow more than one bond between two atoms, so a single
 * bond object must represent both connections.
 */
SKETCHER_API bool contains_two_monomer_linkages(const RDKit::Bond* bond);

} // namespace sketcher
} // namespace schrodinger
