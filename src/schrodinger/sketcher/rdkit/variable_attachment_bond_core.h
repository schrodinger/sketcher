#pragma once

#include <memory>
#include <stdexcept>
#include <unordered_set>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
class RWMol;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{
namespace sketcher
{

/**
 * An exception raised when a variable attachment bond cannot be created because
 * invalid variable attachment atoms have been specified.  (There must be at
 * least two variable attachment atoms and they should be part of the same
 * connected molecule.)
 */
class SKETCHER_API variable_attachment_bond_error : public std::runtime_error
{
  public:
    variable_attachment_bond_error(const std::string& message) :
        std::runtime_error(message){};
};

/**
 * @return Whether the specified bond is a variable attachment bond
 */
SKETCHER_API bool is_variable_attachment_bond(const RDKit::Bond* bond);

/**
 * @return Whether the specified atom is a dummy atom (used to represent the
 * ring center) for a variable attachment bond
 */
SKETCHER_API bool
is_dummy_atom_for_variable_attachment_bond(const RDKit::Atom* atom);

/**
 * @return All variable attachment atoms for the specified variable attachment
 * bond. Will return an empty set if the bond is not a variable attachment bond,
 * or if the variable attachment atom property cannot be parsed.
 */
SKETCHER_API std::unordered_set<const RDKit::Atom*>
get_variable_attachment_atoms(const RDKit::Bond* bond);

/**
 * Calculate coordinates for the specified variable attachment bond, but don't
 * add it to the molecule.
 * @param mol The molecule the bond would be added to. This molecule must have a
 * conformer.
 * @param atoms All atoms that the variable attachment bond should be bound to.
 * This list must contain at least two atoms and the atoms must be part of the
 * same molecule.
 * @return A pair of
 *   - Coordinates for the dummy atom that represents the variable end of the
 *     variable attachment bond
 *   - Coordinates for the carbon atom at the non-variable end of the variable
 *     attachment bond
 */
SKETCHER_API std::pair<RDGeom::Point3D, RDGeom::Point3D>
get_coordinates_for_variable_attachment_bond(
    const RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::Atom*>& atoms);

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
SKETCHER_API std::tuple<RDKit::Atom*, RDKit::Atom*, RDKit::Bond*>
add_variable_attachment_bond_to_mol(
    RDKit::RWMol& mol, const std::unordered_set<const RDKit::Atom*>& atoms);

} // namespace sketcher
} // namespace schrodinger
