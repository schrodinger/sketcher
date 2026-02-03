#pragma once

#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Types of monomers
 */
enum class MonomerType { PEPTIDE, NA_BASE, NA_PHOSPHATE, NA_SUGAR, CHEM };

/**
 * Determine what type of monomer the given atom represents.
 *
 * @throw std::runtime_error if the atom does not represent a monomer
 */
SKETCHER_API MonomerType get_monomer_type(const RDKit::Atom* atom);

/**
 * @return whether the given bond represents two connections between the same
 * monomers, such as neighboring cysteines additionally joined by a disulfide
 * bond. RDKit does not allow more than one bond between two atoms, so a single
 * bond object must represent both connections.
 */
SKETCHER_API bool contains_two_monomer_linkages(const RDKit::Bond* bond);

/**
 * @return a list of {attachment point name, bound atom} for all bound
 * attachment points of the given monomer using "pretty" names (e.g. "N" instead
 * of "R1" for amino acids)
 */
SKETCHER_API std::vector<std::pair<std::string, const RDKit::Atom*>>
get_bound_attachment_point_names_and_atoms(const RDKit::Atom* monomer);

/**
 * @return a list of all unbound attachment points names for the given monomer
 * using "pretty" names (e.g. "N" instead of "R1" for amino acids)
 */
SKETCHER_API std::vector<std::string>
get_available_attachment_point_names(const RDKit::Atom* monomer);

} // namespace sketcher
} // namespace schrodinger
