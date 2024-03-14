#pragma once

#include <string>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

RDKIT_EXTENSIONS_API bool dependency_test(const std::string& smiles);

} // namespace rdkit_extensions
} // namespace schrodinger
