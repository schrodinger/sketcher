#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <memory>
#include <string>

namespace RDKit
{
class RWMol;
}

namespace helm
{
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::RWMol>
helm_to_rdkit(const std::string& helm_string);
}
