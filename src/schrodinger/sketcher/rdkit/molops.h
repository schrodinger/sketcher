#pragma once

#include <string>
#include <unordered_set>
#include <variant>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <GraphMol/ChemReactions/Reaction.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/definitions.h"

using schrodinger::rdkit_extensions::Format;

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

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
 * Create a molecule from a text string.  Atomic coordinates will be
 * automatically generated using coordgen if none are present.
 */
SKETCHER_API boost::shared_ptr<RDKit::RWMol>
text_to_mol(const std::string& text, const Format format = Format::AUTO_DETECT);

/**
 * If the given string represents a molecule, return the specified molecule.  If
 * the string instead represents a reaction, return the specified reaction.
 * @return A tuple of
 *   - If the text represented a molecule, a pointer to the newly created
 *     molecule.  Otherwise nullptr;
 *   - If the text represented a reaction, a pointer to the newly created
 *     reaction.  Otherwise nullptr;
 *   - Whether the text represented a reaction (true) or a molecule (false)
 */
std::variant<boost::shared_ptr<RDKit::RWMol>,
             boost::shared_ptr<RDKit::ChemicalReaction>>
text_to_mol_or_reaction(const std::string& text,
                        const Format format = Format::AUTO_DETECT);

/**
 * Update any RDKit metadata that is required to render the molecule in the
 * Scene.  This function must be called anytime atoms or bonds are changed.
 */
SKETCHER_API void update_molecule_metadata(RDKit::ROMol& mol);

/**
 * @return all atoms and bonds that are connected to the specified atom
 */
SKETCHER_API std::pair<std::unordered_set<const RDKit::Atom*>,
                       std::unordered_set<const RDKit::Bond*>>
get_connected_atoms_and_bonds(const RDKit::Atom* const atom);

} // namespace sketcher
} // namespace schrodinger
