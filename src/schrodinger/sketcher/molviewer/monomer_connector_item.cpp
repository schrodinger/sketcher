#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"

#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/ROMol.h>

#include <QPainter>
#include <QPointF>
#include <QtMath>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

namespace schrodinger
{
namespace sketcher
{

MonomerConnectorItem::MonomerConnectorItem(
    const RDKit::Bond* bond, const AbstractMonomerItem& start_monomer_item,
    const AbstractMonomerItem& end_monomer_item,
    const bool is_secondary_connection, QGraphicsItem* parent) :
    AbstractBondOrConnectorItem(bond, parent),
    m_start_item(start_monomer_item),
    m_end_item(end_monomer_item),
    m_is_secondary_connection(is_secondary_connection)
{

    setZValue(static_cast<qreal>(ZOrder::MONOMER_CONNECTOR));
    updateCachedData();
}

int MonomerConnectorItem::type() const
{
    return Type;
}

bool MonomerConnectorItem::isSecondaryConnection() const
{
    return m_is_secondary_connection;
}

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
static std::tuple<bool, bool, QColor, qreal, Qt::PenStyle>
get_connector_style(const RDKit::Bond* bond, const bool is_secondary_connection)
{
    const auto* start_atom = bond->getBeginAtom();
    auto start_res_name = get_monomer_res_name(start_atom);
    auto start_monomer_type = get_monomer_type(start_atom);

    const auto* end_atom = bond->getEndAtom();
    auto end_res_name = get_monomer_res_name(end_atom);
    auto end_monomer_type = get_monomer_type(end_atom);

    if (start_monomer_type == MonomerType::CHEM ||
        end_monomer_type == MonomerType::CHEM) {
        // chem connector
        return {false, false, CHEM_CONNECTOR_COLOR, CHEM_CONNECTOR_WIDTH,
                Qt::PenStyle::SolidLine};
    } else if (start_monomer_type == MonomerType::PEPTIDE ||
               end_monomer_type == MonomerType::PEPTIDE) {
        std::string attachment_points;
        std::string prop = is_secondary_connection ? CUSTOM_BOND : LINKAGE;
        bond->getPropIfPresent(prop, attachment_points);
        if (start_res_name.ends_with('C') && end_res_name.ends_with('C') &&
            attachment_points == "R3-R3") {
            return {true, true, DISULFIDE_CONNECTOR_COLOR,
                    DISULFIDE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
        }
        bool start_is_branch = false;
        bool end_is_branch = false;
        start_atom->getPropIfPresent(BRANCH_MONOMER, start_is_branch);
        end_atom->getPropIfPresent(BRANCH_MONOMER, end_is_branch);
        if (start_is_branch || end_is_branch) {
            // branching connector
            return {start_is_branch, end_is_branch,
                    AA_BRANCHING_CONNECTOR_COLOR, AA_BRANCHING_CONNECTOR_WIDTH,
                    Qt::PenStyle::SolidLine};
        }
        // standard amino acid linear connector
        return {false, false, AA_LINEAR_CONNECTOR_COLOR,
                AA_LINEAR_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    } else if (start_monomer_type == MonomerType::NA_BASE &&
               end_monomer_type == MonomerType::NA_BASE) {
        // nucleic acid base connector
        return {false, false, NA_BASE_CONNECTOR_COLOR, NA_BASE_CONNECTOR_WIDTH,
                Qt::PenStyle::DotLine};
    } else if (start_monomer_type != MonomerType::NA_BASE &&
               end_monomer_type != MonomerType::NA_BASE) {
        // nucleic acid backbone connector
        return {false, false, NA_BACKBONE_CONNECTOR_COLOR,
                NA_BACKBONE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    } else {
        // nucleic acid backbone to base connector
        return {false, false, NA_BACKBONE_TO_BASE_CONNECTOR_COLOR,
                NA_BACKBONE_TO_BASE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    }
}

/**
 * Return true if coord is above (within a 45 degree cone of) other
 */
static bool is_coord_above_the_other(const QPointF& coord, const QPointF& other)
{
    return coord.y() < other.y() &&
           qFabs(coord.x() - other.x()) < qFabs(coord.y() - other.y());
}

/**
 * Add a diamond shape to the path at the given point
 *
 * @param path the path to add to
 * @param center the center of the diamond
 * @param radius the radius of the diamond
 */
static void add_diamond_arrowhead_to_path(QPainterPath& path,
                                          const QPointF& center,
                                          const qreal radius)
{
    QPolygonF diamond;
    diamond << QPointF(radius, 0) << QPointF(0, -radius) << QPointF(-radius, 0)
            << QPointF(0, radius);
    diamond.translate(center);
    path.addPolygon(diamond);
    path.closeSubpath();
}

/**
 * Update a path so it contains the union of its current value and a diamond
 * shape. If the current contents of path overlap the diamond, then the diamond
 * will be added to an existing subpath instead of creating a new one.
 *
 * @param path the path to add to
 * @param center the center of the diamond
 * @param radius the radius of the diamond
 */
static void or_diamond_arrowhead_to_path(QPainterPath& path,
                                         const QPointF& center,
                                         const qreal radius)
{
    QPainterPath diamond_path;
    add_diamond_arrowhead_to_path(diamond_path, center, radius);
    path |= diamond_path;
}

void MonomerConnectorItem::updateCachedData()
{
    prepareGeometryChange();
    m_arrowhead_path.clear();
    auto [start_has_arrowhead, end_has_arrowhead, connector_color,
          connector_width, connector_pen_style] =
        get_connector_style(m_bond, m_is_secondary_connection);
    m_connector_pen =
        QPen(connector_color, connector_width, connector_pen_style);
    m_arrowhead_pen = QPen(connector_color, connector_width);
    m_arrowhead_pen.setJoinStyle(Qt::PenJoinStyle::MiterJoin);
    m_arrowhead_brush = QBrush(connector_color);

    auto start_qcoords = m_start_item.pos();
    auto end_qcoords = m_end_item.pos();
    setPos(start_qcoords);

    QPointF start_offset;
    if (start_has_arrowhead) {
        auto offset = m_start_item.boundingRect().height() / 2 +
                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
        if (is_coord_above_the_other(start_qcoords, end_qcoords)) {
            offset *= -1;
        }
        start_offset.ry() -= offset;
        add_diamond_arrowhead_to_path(m_arrowhead_path, start_offset,
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    }

    auto end_pos = end_qcoords - start_qcoords;
    if (end_has_arrowhead) {
        auto offset = m_end_item.boundingRect().height() / 2 +
                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
        if (is_coord_above_the_other(end_qcoords, start_qcoords)) {
            offset *= -1;
        }
        end_pos.ry() -= offset;
        add_diamond_arrowhead_to_path(m_arrowhead_path, end_pos,
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    }
    m_connector_line = QLineF(start_offset, end_pos);
    m_midpoint = m_connector_line.center();
    m_selection_highlighting_path = path_around_line(
        m_connector_line, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH);
    m_predictive_highlighting_path = path_around_line(
        m_connector_line, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH);
    if (start_has_arrowhead) {
        or_diamond_arrowhead_to_path(m_selection_highlighting_path,
                                     start_offset,
                                     BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
        or_diamond_arrowhead_to_path(m_predictive_highlighting_path,
                                     start_offset,
                                     BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    }
    if (end_has_arrowhead) {
        or_diamond_arrowhead_to_path(m_selection_highlighting_path, end_pos,
                                     BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
        or_diamond_arrowhead_to_path(m_predictive_highlighting_path, end_pos,
                                     BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    }
    m_shape = QPainterPath(m_selection_highlighting_path);
    m_bounding_rect = m_shape.boundingRect();
}

void MonomerConnectorItem::paint(QPainter* painter,
                                 const QStyleOptionGraphicsItem* option,
                                 QWidget* widget)
{
    painter->save();
    painter->setPen(m_connector_pen);
    painter->drawLine(m_connector_line);

    if (!m_arrowhead_path.isEmpty()) {
        painter->setPen(m_arrowhead_pen);
        painter->setBrush(m_arrowhead_brush);
        painter->drawPath(m_arrowhead_path);
    }

    painter->restore();
}

} // namespace sketcher
} // namespace schrodinger
