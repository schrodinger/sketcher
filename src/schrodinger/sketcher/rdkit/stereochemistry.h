#pragma once

#include <rdkit/GraphMol/GraphMol.h>
#include <QString>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * @return the chiral label for the given atom.
 */
SKETCHER_API QString get_atom_chirality_label(const RDKit::Atom& atom);

/**
 * @return the chiral label for the given bond.
 */
SKETCHER_API QString get_bond_stereo_label(const RDKit::Bond& bond);

/**
 * assign CIP labels to the given molecule
 */
SKETCHER_API void assign_CIP_labels(RDKit::RWMol& mol);

} // namespace sketcher
} // namespace schrodinger