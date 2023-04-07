#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <string>

enum class [[nodiscard]] HelmVersion{V1, V2, Unsupported};

RDKIT_EXTENSIONS_API HelmVersion
get_helm_version(const std::string& helm_string);

RDKIT_EXTENSIONS_API std::string
get_helm2_from_helm(const std::string& helm_string);
