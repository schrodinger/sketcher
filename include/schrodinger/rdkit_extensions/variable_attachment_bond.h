#pragma once

#include <memory>
#include <unordered_set>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace RDKit
{
class Atom;
class Bond;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @return Whether the specified bond is a variable attachment bond
 */
RDKIT_EXTENSIONS_API bool is_variable_attachment_bond(const RDKit::Bond* bond);

/**
 * @return Whether the specified atom is a dummy atom (used to represent the
 * ring center) for a variable attachment bond
 */
RDKIT_EXTENSIONS_API bool
is_dummy_atom_for_variable_attachment_bond(const RDKit::Atom* atom);

/**
 * @return All variable attachment atoms for the specified variable attachment
 * bond. Will return an empty set if the bond is not a variable attachment bond,
 * or if the variable attachment atom property cannot be parsed.
 */
RDKIT_EXTENSIONS_API std::unordered_set<const RDKit::Atom*>
get_variable_attachment_atoms(const RDKit::Bond* bond);

/**
 * Add a variable attachment bond to the specified molecule. The non-variable
 * end of the bond will be bound to a newly-created carbon atom.
 * @param mol The molecule to add the bond to. This molecule must have a
 * conformer.
 * @param atoms All atoms that the variable attachment bond should be bound to.
 * This list must contain at least two atoms and the atoms must be part of the
 * same molecule.
 * @return A tuple of
 *   - The dummy atom that represents the variable end of the variable
 *     attachment bond
 *   - The carbon atom at the non-variable end of the variable attachment bond
 *   - The variable attachment bond
 */
RDKIT_EXTENSIONS_API std::tuple<RDKit::Atom*, RDKit::Atom*, RDKit::Bond*>
add_variable_attachment_bond_to_mol(
    RDKit::RWMol& mol, const std::unordered_set<const RDKit::Atom*>& atoms);

} // namespace rdkit_extensions
} // namespace schrodinger
