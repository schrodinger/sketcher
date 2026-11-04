#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace RDKit
{
class Atom;
class ROMol;
} // namespace RDKit

namespace schrodinger::rdkit_extensions
{

enum class ChainType;

// Read a CXSMILES _R# label, if present.
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<unsigned int>
get_attachment_point_num_from_atom_label(const RDKit::Atom& atom);

// Numbered attachment points may use an atom map, a CXSMILES _R# label, or
// an isotope on a dummy atom, in that order of precedence.
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<unsigned int>
get_attachment_point_num(const RDKit::Atom& atom);

// Return each attachment point number and the heavy element at that site.
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_atomistic_mol(const RDKit::ROMol& mol);

[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_smiles(const std::string& smiles);

// Unknown RNA base N has R1; other missing monomers have R1 and R2.
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_res(const std::string& resname, ChainType chain_type);

// Inspect the inline SMILES or database residue represented by a monomer atom.
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_numbered_attachment_points_for_monomer(const RDKit::Atom* monomer);

// Return the heavy element at one numbered site, if that site exists.
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<std::string>
get_attachment_point_element_for_monomer(const RDKit::Atom* monomer,
                                         unsigned int attachment_point_num);

} // namespace schrodinger::rdkit_extensions
