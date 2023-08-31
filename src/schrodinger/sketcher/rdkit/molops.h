#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include <boost/shared_ptr.hpp>

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
text_to_mol(const std::string& text, Format format = Format::AUTO_DETECT);

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
