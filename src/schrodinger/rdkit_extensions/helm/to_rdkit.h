#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <memory>
#include <string>

namespace RDKit
{
class ROMol;
}

namespace helm
{
RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::ROMol> helm_to_rdkit
    [[nodiscard]] (const std::string& helm_string);
}
