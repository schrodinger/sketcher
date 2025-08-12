#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <string>

namespace RDKit
{
class ROMol;
}

namespace helm
{
RDKIT_EXTENSIONS_API std::string rdkit_to_helm
    [[nodiscard]] (const ::RDKit::ROMol& mol);
} // namespace helm
