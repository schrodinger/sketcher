#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <Qt>

#include "schrodinger/sketcher/definitions.h"

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

enum class Direction { N, S, E, W, NW, NE, SW, SE };

/**
 * Information about an attachment point on a monomer that's bound to another
 * monomer. The direction member variable represents the direction that the bond
 * is drawn in. This is normally the direction of bound_monomer, but may be
 * different if the bond uses an arrowhead (which are typically drawn above or
 * below the monomer).
 */
struct BoundAttachmentPoint {
    std::string name;
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
 * hovers over the monomer.
 */
struct UnboundAttachmentPoint {
    std::string name;
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
 * distance between the center of the monomer and the center of the arrowhead.
 * @param monomer_item the graphics item for the monomer
 * @param bound_coords the Scene coordinates for the other monomer involved in
 * the bond
 */
SKETCHER_API qreal get_monomer_arrowhead_offset(
    const QGraphicsItem& monomer_item, const QPointF& bound_coords);

/**
 * @return a list of all unbound attachment points names for the given monomer
 * using "pretty" names (e.g. "N" instead of "R1" for amino acids)
 */

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

} // namespace sketcher
} // namespace schrodinger
