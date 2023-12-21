#pragma once
#include <unordered_set>

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a set of the atoms in the smaller of the two subset that the given
 * bond divides the molecule into. If the bond is part of a ring, an empty set
 * is returned.
 * @param mol The molecule to be split
 * @param bond The bond that divides the molecule into two subsets
 */
std::unordered_set<const RDKit::Atom*>
get_smaller_substituent_atoms(const RDKit::ROMol& mol, const RDKit::Bond& bond);

} // namespace sketcher
} // namespace schrodinger
