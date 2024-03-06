#pragma once

#include <memory>

#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * @return a dummy atom, properly configured so that it will be represented
 * correctly in SMARTS and MDL MOL file output
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::shared_ptr<RDKit::Atom>
create_dummy_atom();

} // namespace rdkit_extensions
} // namespace schrodinger
