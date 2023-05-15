#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <string>
#include <memory> // shared_ptr

namespace RDKit
{
class ROMol;
} // namespace RDKit

RDKIT_EXTENSIONS_API std::shared_ptr<::RDKit::ROMol>
helm_to_rdkit(const std::string& input_helm);
