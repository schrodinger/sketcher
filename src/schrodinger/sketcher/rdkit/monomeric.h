#pragma once

#include <concepts>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include <fmt/format.h>
#include <Qt>

#include "schrodinger/rdkit_extensions/monomer_directions.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/definitions.h"

class QGraphicsItem;
class QPointF;

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{
namespace sketcher
{

using rdkit_extensions::Direction;

enum class MonomerType { PEPTIDE, NA_BASE, NA_PHOSPHATE, NA_SUGAR, CHEM };

enum class ConnectorType {
    CHEM,
    PEPTIDE_LINEAR,
    PEPTIDE_BRANCHING,
    PEPTIDE_DISULFIDE,
    PEPTIDE_SIDE_CHAIN,
    H_BOND,
    NA_BACKBONE,
    NA_BACKBONE_TO_BASE
};

namespace PeptideAP
{
enum { N = 1, C = 2, X_OR_S = 3 };
}

namespace NASugarAP
{
enum { FIVE_PRIME = 1, THREE_PRIME = 2, ONE_PRIME = 3 };
}

namespace NAPhosphateAP
{
enum { TO_PREV_SUGAR = 1, TO_NEXT_SUGAR = 2 };
}

constexpr int NA_BASE_AP_N1_9 = 1;
const std::string H_BOND_AP_MODEL_NAME = "pair";
const std::string H_BOND_DISPLAY_NAME = "H-bond";

const std::string PEPTIDE_R3_NAME_S = "S";
const std::string PEPTIDE_R3_NAME_X = "X";

/**
 * Validate the monomers in a monomeric molecule.
 * @throw std::runtime_error if a monomer is missing from the database (other
 * than peptide X and nucleic acid N, which represent unknown monomers) or its
 * inline SMILES cannot be parsed.
 */
SKETCHER_API void validate_monomers(const RDKit::ROMol& mol);

/**
 * Convert any of the above attachment point enums (or NA_BASE_AP_N1_9) to the
 * equivalent model name, which is simply "R" followed by the attachment point
 * number.
 */
SKETCHER_API std::string ap_model_name_for(int ap_num);

/**
 * @return a list of all attachment points in required_attachment_points that
 * are not present in atomistic_mol
 */
SKETCHER_API std::vector<int> get_missing_required_attachment_points(
    const RDKit::ROMol& atomistic_mol,
    std::vector<int> required_attachment_points);

/**
 * Return the numbered attachment points in the given atomistic molecule (which
 * contains the contents of a monomer). Each attachment point is described using
 * a pair of the attachment point number and the symbol of the heavy atom at
 * that site. Note that this function assumes that the attachment points in the
 * molecule are sane; it does not protect against, e.g., duplicated attachment
 * points or attachment points on unbound dummy atoms.
 *
 * @throws std::invalid_argument if smiles is not valid extended SMILES
 */
SKETCHER_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_atomistic_mol(const RDKit::ROMol& mol);

/**
 * Return the numbered attachment points in a monomer SMILES string. Each
 * attachment point is described using a pair of the attachment point number and
 * the symbol of the heavy atom at that site. Note that this function assumes
 * that the SMILES string is valid and sane; it does not protect against, e.g.,
 * duplicated attachment points or attachment points on unbound dummy atoms.
 *
 * @param smiles the SMILES string, with attachment points indicated using
 * atom-map numbers, isotope-numbered dummy atoms, or CXSMILES atom labels such
 * as "_R1". For example, any of the following are acceptable descriptions of
 * alanine:
 *   - C[C@H](N[H:1])C(=O)[OH:2]
 *   - [1*]N[C@@H](C)C(=O)O[2*]
 *   - *N[C@@H](C)C(=O)O* |$_R1;;;;;;;_R2$|.
 *
 * @throws std::invalid_argument if smiles is not valid extended SMILES
 */
SKETCHER_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_smiles(const std::string& smiles);

/**
 * Normalize numbered attachment points in a monomer SMILES string so each is
 * represented by a dummy atom with a CXSMILES atom label such as "_R1".
 * Attachment points may initially be represented by atom-map numbers,
 * isotope-numbered dummy atoms, or CXSMILES atom labels. Atom-mapped leaving
 * atoms are replaced by labeled dummy atoms. A heavy atom marked only with a
 * CXSMILES label receives a new bonded dummy atom.
 *
 * @throws std::invalid_argument if smiles is not valid extended SMILES
 */
SKETCHER_API std::string
normalize_smiles_attachment_points(const std::string& smiles);

/**
 * Return the numbered attachment points for the given monomer. Each attachment
 * point is described using a pair of the attachment point number and the symbol
 * of the heavy atom at that site.
 *
 * If the monomer is missing from the database, unknown nucleic acid bases (N)
 * have R1, while unknown peptides (X) and all other missing monomers have R1
 * and R2. These fallback attachment points have empty element symbols.
 */
SKETCHER_API std::vector<std::pair<int, std::string>>
get_attachment_points_for_res(const std::string& resname,
                              const rdkit_extensions::ChainType chain_type);

/**
 * Information about an attachment point on a monomer that's bound to another
 * monomer. The direction member variable represents the direction that the bond
 * is drawn in. This is normally the direction of bound_monomer, but may be
 * different if the bond uses an arrowhead (which are typically drawn above or
 * below the monomer).
 */
struct BoundAttachmentPoint {
    /// The attachment point name stored in the model. Typically "R" followed by
    /// a positive, non-zero integer
    std::string model_name;
    /// The attachment point name that is displayed in the Sketcher. E.g. "C"
    /// instead of "R2" for a peptide monomer.  Note that, for an attachment
    /// point with a custom name, this will be identical to model_name.
    std::string display_name;
    /// The attachment point number. For attachment points with model names of
    /// an "R" followed by a number, this will be the number that appears after
    /// the "R" (e.g. 3 for "R3"). For attachment points with a custom name,
    /// this will be ATTACHMENT_POINT_WITH_CUSTOM_NAME.
    int num;
    const RDKit::Atom* bound_monomer;
    bool is_secondary_connection;
    Direction direction;

    bool operator==(const BoundAttachmentPoint&) const = default;
};

/**
 * Information about an attachment point on a monomer that's *not* bound to
 * another monomer (i.e. available for bonding). The direction member variable
 * represents the direction we should draw the connection "nubbin" when the user
 * hovers over the monomer. See BoundAttachmentPoint for documentation of member
 * variables.
 */
struct UnboundAttachmentPoint {
    std::string model_name;
    std::string display_name;
    int num;
    Direction direction;

    bool operator==(const UnboundAttachmentPoint&) const = default;
};

/**
 * For both BoundAttachmentPoints and UnboundAttachmentPoints, num will be
 * ATTACHMENT_POINT_WITH_CUSTOM_NAME if the attachment uses a name that doesn't
 * follow the standard R# pattern (e.g. the "pair" attachment point on nucleic
 * acid bases).
 */
const int ATTACHMENT_POINT_WITH_CUSTOM_NAME = -1;

/**
 * Determine what type of monomer the given atom represents.
 *
 * @throw std::runtime_error if the atom does not represent a monomer
 */
SKETCHER_API MonomerType get_monomer_type(const RDKit::Atom* atom);

/**
 * @return the type of nucleic acid monomer that the given residue name
 * represents
 */
SKETCHER_API MonomerType
get_na_monomer_type_from_res_name(const std::string_view res_name);

/**
 * Determine the text to use for the name of the given monomer
 */
SKETCHER_API std::string get_monomer_res_name(const RDKit::Atom* const monomer);

/**
 * @return true if the given NA base atom is part of a DNA strand, i.e. its
 * bound sugar's residue name is "dR". Returns false for RNA bases (sugar "R")
 * and for bases with no bound sugar or a non-standard sugar.
 *
 * @throw std::runtime_error if the atom does not represent an NA_BASE monomer.
 */
SKETCHER_API bool is_dna_base(const RDKit::Atom* const base);

/**
 * @return the Watson-Crick DNA complement symbol for the given nucleic acid
 * base residue name, or std::nullopt if the symbol has no standard complement.
 * @param base_symbol the residue name of the base (e.g. "A", "G")
 */
SKETCHER_API std::optional<std::string>
get_dna_complement_base_symbol(const std::string_view base_symbol);

/**
 * @return the Watson-Crick RNA complement symbol for the given nucleic acid
 * base residue name, or std::nullopt if the symbol has no standard complement.
 * @param base_symbol the residue name of the base (e.g. "A", "G")
 */
SKETCHER_API std::optional<std::string>
get_rna_complement_base_symbol(const std::string_view base_symbol);

/**
 * @return true if the given nucleic acid base residue name has a standard
 * Watson-Crick complement. Existence is independent of the DNA/RNA target, so
 * this is suitable for gating UI before the target strand type is known.
 */
SKETCHER_API bool na_base_has_complement(const std::string_view base_symbol);

/**
 * The data needed to build one nucleotide (sugar, base, and phosphate) of a
 * complementary strand
 */
struct ComplementNucleotide {
    // the index of the original base that this nucleotide pairs with
    size_t original_base_idx;
    std::string sugar_symbol;
    std::string base_symbol;
};

/**
 * Determine the complementary chains needed to pair with the given nucleic acid
 * bases. Bases are grouped by polymer, and each polymer's bases are split into
 * runs of neighboring nucleotides, each of which gets its own complementary
 * chain. Bases without a Watson-Crick complement are skipped (which also splits
 * the run), as are atoms that aren't nucleic acid bases.
 *
 * @param bases the nucleic acid bases to complement
 * @return the complementary chains in the order that they should be created.
 * Each chain lists its nucleotides in the residue order of the original bases
 * that they pair with.
 */
SKETCHER_API std::vector<std::vector<ComplementNucleotide>>
get_complement_chains(const std::unordered_set<const RDKit::Atom*>& bases);

/**
 * Determine the direction from the original bases toward their complements. The
 * complement sits on the side of the base away from the base's own sugar, so
 * deriving the direction from the geometry lets the complement follow a rotated
 * or moved strand.
 *
 * @param mol the monomer molecule containing the original bases
 * @param complement_nucleotides the nucleotides of one complementary chain
 * @return a unit vector toward the complement, taken from the first original
 * base with a locatable sugar, or (0, -1) if there is no such base
 */
SKETCHER_API RDGeom::Point3D get_complement_pairing_direction(
    const RDKit::ROMol& mol,
    const std::vector<ComplementNucleotide>& complement_nucleotides);

/**
 * @return whether the given bond represents two connections between the same
 * monomers, such as neighboring cysteines additionally joined by a disulfide
 * bond. RDKit does not allow more than one bond between two atoms, so a single
 * bond object must represent both connections.
 */
SKETCHER_API bool contains_two_monomer_linkages(const RDKit::Bond* bond);

/**
 * Determine what type of monomeric connection the given bond represents.
 * Note that a connector between a CHEM monomer and any other monomer will be
 * categorized as a CHEM connector, and a connector between a PEPTIDE monomer
 * and an RNA monomer will be categorized as a PEPTIDE connector.
 *
 * @param bond the bond representing a monomer connector
 * @param is_secondary_connection whether to categorize the primary or secondary
 * connection of this bond. This is only relevant for bonds that
 * contains_two_monomer_linkages() returns true.
 */
SKETCHER_API ConnectorType
get_connector_type(const RDKit::Bond* bond, const bool is_secondary_connection);

/**
 * Determine whether diamond arrowheads should be drawn at the start and end of
 * the given connector.
 * @param bond the bond representing a monomer connector
 * @param is_secondary_connection whether we are drawing the primary or
 * secondary connection of this bond. This is only relevant for bonds where
 * contains_two_monomer_linkages() returns true.
 * @return A pair of
 *   - should an arrowhead be drawn at the start of the bond
 *   - should an arrowhead be drawn at the end of the bond
 */
std::pair<bool, bool>
does_connector_have_arrowheads(const RDKit::Bond* bond,
                               const bool is_secondary_connection);

/**
 * @overload
 * @param bond the bond representing a monomer connector
 * @param connector_type the type of monomer connection represented by bond
 */
std::pair<bool, bool>
does_connector_have_arrowheads(const RDKit::Bond* bond,
                               const ConnectorType connector_type);

/**
 * For a monomer connector being drawn with a diamond arrowhead, determine the
 * offset from the center of the monomer to the center of the arrowhead. The
 * arrowhead is placed outside the nearest cardinal side facing the bound
 * monomer. If that side is occupied, up to three other forward-facing sides
 * and corners are considered, preferring sides over corners. A side or corner
 * is occupied when another connection is within 45 degrees of its direction.
 * Independently placed endpoints use the candidate crossing the fewest other
 * bonds, with normal direction priority breaking ties. Same-chain endpoints
 * are coordinated when necessary so their connector does not cross the chain
 * between them.
 * @param monomer_item the graphics item for the monomer
 * @param bound_coords the Scene coordinates for the other monomer involved in
 * the bond
 * @param monomer the monomer where the arrowhead is being placed
 * @param bound_monomer the other monomer involved in the connection
 * @param is_secondary_connection whether this is the secondary connection of
 * a bond that represents two monomer connections
 */
SKETCHER_API QPointF get_monomer_arrowhead_offset(
    const QGraphicsItem& monomer_item, const QPointF& bound_coords,
    const RDKit::Atom* monomer, const RDKit::Atom* bound_monomer,
    bool is_secondary_connection);

/**
 * Low-level overload that uses an explicitly provided set of occupied sides
 * and corners. This is useful when the caller has already determined
 * connection occupancy.
 * @param monomer_item the graphics item for the monomer
 * @param bound_coords the Scene coordinates for the other monomer
 * @param occupied_directions sides and corners already used by other
 * connections
 */
SKETCHER_API QPointF get_monomer_arrowhead_offset(
    const QGraphicsItem& monomer_item, const QPointF& bound_coords,
    const std::unordered_set<Direction>& occupied_directions);

/**
 * @return whether or not the specified peptide monomer has a side-chain
 * attachment point
 * @param res_name_or_smiles The residue name if is_smiles is false or the
 * SMILES string if is_smiles if true.
 * @param is_smiles Whether res_name_or_smiles represents the residue name (for
 * residues that appear in the monomer database) or the SMILES string (for
 * SMILES monomers)
 */
SKETCHER_API bool peptide_has_ap3(const std::string& res_name_or_smiles,
                                  const bool is_smiles);

/**
 * Determine all bound and unbound monomeric attachment points for the given
 * monomer
 * @return A pair of
 *   - A list of all bound attachment points containing "pretty" names (e.g. "N"
 *     instead of "R1" for amino acids). The direction represents the direction
 *     that the bond is drawn in.
 *   - A list of all unbound attachment points containing "pretty" names (e.g.
 *     "N" instead of "R1" for amino acids). The direction represents the
 *     direction that the attachment point indicator should be drawn.
 */
SKETCHER_API std::pair<std::vector<BoundAttachmentPoint>,
                       std::vector<UnboundAttachmentPoint>>
get_attachment_points_for_monomer(const RDKit::Atom* monomer);

/**
 * Return the attachment point name for the specified monomeric connection
 * @param monomer the monomer to determine the attachment point for
 * @param connector the monomeric connection to get the attachment point for
 * @param is_secondary_connection if this name is for the secondary connection
 * of the bond
 * @return the "pretty" attachment point name (e.g. "N" instead of "R1" for
 * amino acids)
 */
SKETCHER_API std::string
get_attachment_point_name_for_connection(const RDKit::Atom* monomer,
                                         const RDKit::Bond* connector,
                                         const bool is_secondary_connection);

/**
 * @return the attachment point number represented by the given name, e.g. 3 for
 * "R3". If the attachment point name doesn't follow the standard "R" followed
 * by a number format, then ATTACHMENT_POINT_WITH_CUSTOM_NAME will be returned.
 */
SKETCHER_API int ap_name_to_num(const std::string_view attachment_point_name);

/**
 * Take all monomers in the merge_from chain/polymer, and add them to the
 * merge_to chain/polymer.
 */
SKETCHER_API void merge_chains(RDKit::ROMol& mol,
                               const std::string_view merge_from,
                               const std::string& merge_to);

/**
 * Return the numeric component of the chain name as an integer, e.g. 3 for
 * "PEPTIDE3". If the chain name cannot be parsed, -1 will be returned.
 * @param chain_name the chain name to parse
 * @param chain_type the type of chain_name
 */
SKETCHER_API int get_chain_num(const std::string_view chain_name,
                               const rdkit_extensions::ChainType chain_type);

/**
 * @return the lowest numbered chain name (of the specified type) that doesn't
 * already exist in the molecule. For example, if a molecule already has
 * "PEPTIDE1" and "PEPTIDE2" chains, "PEPTIDE3" will be returned when
 * ChainType::PEPTIDE is passed in. Alternatively, "RNA1" would be returned if
 * ChainType::RNA was passed in.
 */
SKETCHER_API std::string
get_first_available_chain_name(const RDKit::ROMol& mol,
                               const rdkit_extensions::ChainType chain_type);

/**
 * Combine the two attachment point names to form a standardized linkage string.
 * For example, attachment points "R2" and "R1" would form the linkage string
 * "R2-R1". Note that, in a standardized linkage string, higher numbered
 * attachment points are listed before lower numbered one.
 *
 * Also note that the "pair" attachment point implies hydrogen bonding, while
 * the "R<#>" attachment points imply covalent or disulfide bonding, so a "pair"
 * attachment point can't be connected to a numbered attachment point.  As a
 * result, a linkage of "pair-pair" will be returned if either attachment point
 * is "pair".
 * @return a pair of
 *   - the standardized linkage string
 *   - whether the attachment point order was flipped in order to standardize
 *     the linkage string
 */
SKETCHER_API std::pair<std::string, bool>
build_linkage_string(const std::string_view ap_name_one,
                     const std::string_view ap_name_two);

/**
 * Determine whether the described linkage is a standard bond (i.e. does not
 * need the CUSTOM_BOND property) or a custom bond (which requires the
 * CUSTOM_BOND) property
 * @note the linkage string must be standardized, but the order of the monomers
 * does not need to match the order of the linkage string. In a standardized
 * linkage string, higher numbered attachment points are listed before lower
 * numbered one, and numbered attachment points are listed before attachment
 * points with custom names.
 */
SKETCHER_API bool
get_is_custom_bond(const std::string_view res_name_one,
                   const rdkit_extensions::ChainType chain_type_one,
                   const std::string_view res_name_two,
                   const rdkit_extensions::ChainType chain_type_two,
                   const std::string_view linkage);

/// @overload
SKETCHER_API bool get_is_custom_bond(const RDKit::Atom* const monomer_one,
                                     const RDKit::Atom* const monomer_two,
                                     const std::string_view linkage);

/// @overload
SKETCHER_API bool
get_is_custom_bond(const std::string_view res_name_one,
                   const rdkit_extensions::ChainType chain_type_one,
                   const RDKit::Atom* const monomer_two,
                   const std::string_view linkage);
/**
 * Add a connection between the two given monomers and their respective
 * attachment points.
 * @return The index of the bond containing the newly added connection.
 */
SKETCHER_API unsigned int
add_monomer_connection(RDKit::RWMol& mol, const unsigned int monomer_one_idx,
                       const std::string_view ap_one,
                       const unsigned int monomer_two_idx,
                       const std::string_view ap_two);

/// @overload
SKETCHER_API unsigned int add_monomer_connection(RDKit::RWMol& mol,
                                                 RDKit::Atom* const monomer_one,
                                                 const std::string_view ap_one,
                                                 RDKit::Atom* const monomer_two,
                                                 const std::string_view ap_two);

/**
 * Change chain names in mol_to_change to ensure that all of the names are
 * different from those in reference_mol
 */
SKETCHER_API void
ensure_distinct_chain_names(RDKit::ROMol& mol_to_change,
                            const RDKit::ROMol& reference_mol);

/**
 * @return whether the given bond represents a monomeric hydrogen bond
 */
SKETCHER_API bool is_hydrogen_bond(const RDKit::Bond* const bond);

} // namespace sketcher
} // namespace schrodinger
