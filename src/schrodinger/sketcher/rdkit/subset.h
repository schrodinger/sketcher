#pragma once
#include <unordered_set>

#include "schrodinger/sketcher/definitions.h"

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
 * @return whether the given atoms and bonds form a contiguous region in the
 * molecule. i.e. whether the atoms are connected to each other through the
 * bonds
 */
SKETCHER_API
bool is_contiguous_region(std::unordered_set<const RDKit::Atom*> atoms,
                          std::unordered_set<const RDKit::Bond*> bonds);

/**
 * @return all atoms and bonds that are connected to the specified atom
 */
SKETCHER_API std::pair<std::unordered_set<const RDKit::Atom*>,
                       std::unordered_set<const RDKit::Bond*>>
get_connected_atoms_and_bonds(const RDKit::Atom* const atom);

/**
 * @return whether the atoms are part of the same fragment or not
 */
SKETCHER_API bool
in_same_fragment(const std::unordered_set<const RDKit::Atom*>& atoms);

/**
 * @return a set of the atoms in the smaller of the two subset that the given
 * bond divides the molecule into. If the bond is part of a ring, an empty set
 * is returned.
 * @param mol The molecule to be split
 * @param bond The bond that divides the molecule into two subsets
 */
SKETCHER_API std::unordered_set<const RDKit::Atom*>
get_smaller_substituent_atoms(const RDKit::ROMol& mol, const RDKit::Bond& bond);

} // namespace sketcher
} // namespace schrodinger
