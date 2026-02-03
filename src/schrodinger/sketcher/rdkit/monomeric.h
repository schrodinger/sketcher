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

/**
 * Types of monomers
 */
enum class MonomerType { PEPTIDE, NA_BASE, NA_PHOSPHATE, NA_SUGAR, CHEM };

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
 * Determine how the given monomer connector should be drawn, which is based on
 * the types of monomers that it connects. Note that a connector between a CHEM
 * monomer and any other monomer will be styled as a CHEM connector, and a
 * connector between a PEPTIDE monomer and an RNA monomer will be styled as a
 * PEPTIDE connector.
 *
 * @param bond the bond representing a monomer connector
 * @param secondary_connection whether we are drawing the primary or secondary
 * connection of this bond
 * @return A tuple of
 *   - whether the start of the connector should have a diamond arrowhead
 *   - whether the end of the connector should have a diamond arrowhead
 *   - the color for the connector
 *   - the width of the connector line
 *   - the pen style of the connector line
 */
SKETCHER_API std::tuple<bool, bool, QColor, qreal, Qt::PenStyle>
get_connector_style(const RDKit::Bond* bond,
                    const bool is_secondary_connection);

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
