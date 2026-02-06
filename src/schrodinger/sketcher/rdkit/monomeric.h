#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <Qt>

#include "schrodinger/sketcher/definitions.h"

class QColor;
class QGraphicsItem;
class QPointF;

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

enum class MonomerType { PEPTIDE, NA_BASE, NA_PHOSPHATE, NA_SUGAR, CHEM };

enum class ConnectorType {
    CHEM,
    PEPTIDE_LINEAR,
    PEPTIDE_BRANCHING,
    PEPTIDE_DISULFIDE,
    PEPTIDE_SIDE_CHAIN,
    NA_BASE,
    NA_BACKBONE,
    NA_BACKBONE_TO_BASE
};

/**
 * Determine what type of monomer the given atom represents.
 *
 * @throw std::runtime_error if the atom does not represent a monomer
 */
SKETCHER_API MonomerType get_monomer_type(const RDKit::Atom* atom);

/**
 * Determine the text to use for the name of the given monomer
 */
SKETCHER_API std::string get_monomer_res_name(const RDKit::Atom* const monomer);

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
 * @param secondary_connection whether to categorize the primary or secondary
 * connection of this bond. This is only relevant for bonds that
 * contains_two_monomer_linkages() returns true.
 */
SKETCHER_API ConnectorType
get_connector_type(const RDKit::Bond* bond, const bool is_secondary_connection);

/**
 * Determine whether diamond arrowheads should be drawn at the start and end of
 * the given connector.
 * @param bond the bond representing a monomer connector
 * @param secondary_connection whether we are drawing the primary or secondary
 * connection of this bond. This is only relevant for bonds where
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
 * distance between the center of the monomer and the center of the arrowhead.
 * @param monomer_item the graphics item for the monomer
 * @param bound_coords the Scene coordinates for the other monomer involved in
 * the bond
 */
SKETCHER_API qreal get_monomer_arrowhead_offset(
    const QGraphicsItem& monomer_item, const QPointF& bound_coords);

/**
 * @return a list of {attachment point name, bound atom, is this the secondary
 * connection of the RDKit::Bond} for all bound attachment points of the given
 * monomer using "pretty" names (e.g. "N" instead of "R1" for amino acids)
 */
SKETCHER_API std::vector<std::tuple<std::string, const RDKit::Atom*, bool>>
get_bound_attachment_point_names_and_atoms(const RDKit::Atom* monomer);

/**
 * @return a list of all unbound attachment points names for the given monomer
 * using "pretty" names (e.g. "N" instead of "R1" for amino acids)
 */
SKETCHER_API std::vector<std::string>
get_available_attachment_point_names(const RDKit::Atom* monomer);

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
get_attachment_point_name_for_atom(const RDKit::Atom* monomer,
                                   const RDKit::Bond* connector,
                                   const bool is_secondary_connection);

} // namespace sketcher
} // namespace schrodinger
