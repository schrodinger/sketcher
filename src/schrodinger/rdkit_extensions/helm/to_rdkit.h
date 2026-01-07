#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <memory>
#include <string>
#include <string_view>

namespace RDKit
{
class RWMol;
}

namespace helm
{
/**
 * @brief Converts a HELM string to an RDKit molecule.
 *
 * Parses a HELM notation string and builds a corresponding RDKit `RWMol`.
 *
 * @param helm_string  Input HELM string (e.g. "BLOB1{Bead}$$$$V2.0").
 * @param do_throw  If true, throw on parse errors; otherwise log and return
 *        nullptr (default: true).
 *
 * @return Parsed RDKit molecule, or nullptr on failure when
 *         `do_throw == false`.
 * @throws std::invalid_argument on invalid HELM input when `do_throw == true`.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<::RDKit::RWMol>
helm_to_rdkit(const std::string& helm_string, bool do_throw = true);

} // namespace helm
