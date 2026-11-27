#pragma once

/**
 * Numbered attachment points on HELM monomers. These helpers inspect the
 * atomistic SMILES definition of a monomer to identify its available sites;
 * the RDKit atoms passed to the monomer helpers still represent whole monomers.
 */

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

/**
 * Read a positive attachment point number from a CXSMILES atom label such as
 * "_R3". Return std::nullopt when no positive number can be read.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<unsigned int>
get_attachment_point_num_from_atom_label(const RDKit::Atom& atom);

/**
 * Get an atom's numbered attachment point, if any. Recognize a positive atom
 * map number first, then a CXSMILES "_R#" label, then an isotope number on a
 * dummy atom. Return std::nullopt if none is present.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<unsigned int>
get_attachment_point_num(const RDKit::Atom& atom);

/**
 * Return the numbered attachment points in an atomistic monomer molecule.
 * Each pair contains the attachment point number and the symbol of its heavy
 * atom. This assumes valid attachment points; duplicate numbers and unbound
 * dummy atoms are not checked.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_atomistic_mol(const RDKit::ROMol& mol);

/**
 * Return numbered attachment points and their heavy atom symbols from extended
 * SMILES. Attachment points may be atom-map numbers, isotope-numbered dummy
 * atoms, or CXSMILES labels such as "_R1". For example, all three forms below
 * describe alanine:
 *   C[C@H](N[H:1])C(=O)[OH:2]
 *   [1*]N[C@@H](C)C(=O)O[2*]
 *   *N[C@@H](C)C(=O)O* |$_R1;;;;;;;_R2$|
 *
 * As with the molecule overload, duplicate or unbound attachment points are
 * not validated.
 * @throws std::invalid_argument if smiles is not valid extended SMILES.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_smiles(const std::string& smiles);

/**
 * Return numbered attachment points and heavy atom symbols for a database
 * monomer. If it is missing, unknown RNA base N has R1; all other missing
 * monomers have R1 and R2. Fallback points have empty element symbols.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_res(const std::string& resname, ChainType chain_type);

/**
 * Return numbered attachment points for a whole-monomer atom. Inline SMILES
 * monomers are parsed from their atom label; named monomers are looked up using
 * their PDB residue name and HELM chain type.
 * @throws std::invalid_argument if a named monomer has no PDB residue info.
 * @throws std::out_of_range if its polymer ID has no recognized HELM type.
 * An invalid inline SMILES string also propagates the parser's exception.
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<std::pair<int, std::string>>
get_numbered_attachment_points_for_monomer(const RDKit::Atom* monomer);

/**
 * Return the heavy atom symbol at one numbered attachment point, or
 * std::nullopt if the point is absent. Propagates the exceptions from
 * get_numbered_attachment_points_for_monomer().
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::optional<std::string>
get_attachment_point_element_for_monomer(const RDKit::Atom* monomer,
                                         unsigned int attachment_point_num);

} // namespace schrodinger::rdkit_extensions
