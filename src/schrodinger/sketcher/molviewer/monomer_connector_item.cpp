#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"

#include <cstdlib>

#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/ROMol.h>

#include <QPainter>
#include <QPointF>
#include <QTransform>
#include <QtMath>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

// clang-format off
std::unordered_map<ConnectorType, std::tuple<QColor, QColor, qreal, Qt::PenStyle>>
    STYLE_FOR_CONNECTOR_TYPE = {
        {ConnectorType::CHEM,
         {CHEM_CONNECTOR_COLOR,                CHEM_CONNECTOR_COLOR_DARK_BG,
          CHEM_CONNECTOR_WIDTH,                Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_LINEAR,
         {AA_LINEAR_CONNECTOR_COLOR,           AA_LINEAR_CONNECTOR_COLOR_DARK_BG,
          AA_LINEAR_CONNECTOR_WIDTH,           Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_BRANCHING,
         {AA_BRANCHING_CONNECTOR_COLOR,        AA_BRANCHING_CONNECTOR_COLOR_DARK_BG,
          AA_BRANCHING_CONNECTOR_WIDTH,        Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_DISULFIDE,
         {DISULFIDE_CONNECTOR_COLOR,           DISULFIDE_CONNECTOR_COLOR_DARK_BG,
          DISULFIDE_CONNECTOR_WIDTH,           Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_SIDE_CHAIN,
         {AA_LINEAR_CONNECTOR_COLOR,           AA_LINEAR_CONNECTOR_COLOR_DARK_BG,
          AA_LINEAR_CONNECTOR_WIDTH,           Qt::PenStyle::SolidLine}},
        {ConnectorType::H_BOND,
         {H_BOND_CONNECTOR_COLOR,             H_BOND_CONNECTOR_COLOR_DARK_BG,
          H_BOND_CONNECTOR_WIDTH,             Qt::PenStyle::DotLine}},
        {ConnectorType::NA_BACKBONE,
         {NA_BACKBONE_CONNECTOR_COLOR,         NA_BACKBONE_CONNECTOR_COLOR_DARK_BG,
          NA_BACKBONE_CONNECTOR_WIDTH,         Qt::PenStyle::SolidLine}},
        {ConnectorType::NA_BACKBONE_TO_BASE,
         {NA_BACKBONE_TO_BASE_CONNECTOR_COLOR, NA_BACKBONE_TO_BASE_CONNECTOR_COLOR_DARK_BG,
          NA_BACKBONE_TO_BASE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine}}};
// clang-format on

} // namespace

MonomerConnectorItem::MonomerConnectorItem(
    const RDKit::Bond* bond, const AbstractMonomerItem& start_monomer_item,
    const AbstractMonomerItem& end_monomer_item,
    const bool is_secondary_connection, const bool is_dark_mode,
    const int connector_lane, QGraphicsItem* parent) :
    AbstractBondOrConnectorItem(bond, parent),
    m_is_dark_mode(is_dark_mode),
    m_start_item(start_monomer_item),
    m_end_item(end_monomer_item),
    m_is_secondary_connection(is_secondary_connection),
    m_connector_lane(connector_lane)
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

int MonomerConnectorItem::getConnectorLane() const
{
    return m_connector_lane;
}

/**
 * Add a diamond shape to the path at the given point
 *
 * @param path the path to add to
 * @param center the center of the diamond
 * @param radius the radius of the diamond
 * @param rotate_for_corner whether to rotate the diamond 45 degrees
 */
static void add_diamond_arrowhead_to_path(QPainterPath& path,
                                          const QPointF& center,
                                          const qreal radius,
                                          const bool rotate_for_corner)
{
    QPolygonF diamond;
    diamond << QPointF(radius, 0) << QPointF(0, -radius) << QPointF(-radius, 0)
            << QPointF(0, radius);
    if (rotate_for_corner) {
        QTransform transform;
        transform.rotate(45.0);
        diamond = transform.map(diamond);
    }
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
 * @param rotate_for_corner whether to rotate the diamond 45 degrees
 */
static void or_diamond_arrowhead_to_path(QPainterPath& path,
                                         const QPointF& center,
                                         const qreal radius,
                                         const bool rotate_for_corner)
{
    QPainterPath diamond_path;
    add_diamond_arrowhead_to_path(diamond_path, center, radius,
                                  rotate_for_corner);
    path |= diamond_path;
}

static bool offset_is_at_corner(const QPointF& offset)
{
    return !qFuzzyIsNull(offset.x()) && !qFuzzyIsNull(offset.y());
}

void MonomerConnectorItem::updateCachedData()
{
    prepareGeometryChange();
    m_arrowhead_path.clear();
    m_start_join_line = QLineF();
    m_end_join_line = QLineF();
    auto connector_type = get_connector_type(m_bond, m_is_secondary_connection);
    auto [start_has_arrowhead, end_has_arrowhead] =
        does_connector_have_arrowheads(m_bond, connector_type);
    auto [connector_color, connector_color_dark_bg, connector_width,
          connector_pen_style] = STYLE_FOR_CONNECTOR_TYPE.at(connector_type);
    m_connector_color = connector_color;
    m_connector_color_dark_bg = connector_color_dark_bg;
    auto color = getConnectorColor();
    m_connector_pen = QPen(color, connector_width, connector_pen_style);
    m_arrowhead_pen = QPen(color, connector_width);
    m_arrowhead_pen.setJoinStyle(Qt::PenJoinStyle::MiterJoin);
    m_arrowhead_brush = QBrush(color);

    auto start_qcoords = m_start_item.pos();
    auto end_qcoords = m_end_item.pos();
    setPos(start_qcoords);

    // A nonzero lane identifies a custom connection in a linear peptide. Its
    // endpoint anchors are fixed above or below the residue boxes. Unrouted
    // connectors (lane zero) keep the normal direction-aware geometry.
    QPointF start_offset{};
    QPointF end_offset{};
    if (m_connector_lane == 0) {
        if (start_has_arrowhead) {
            start_offset = get_monomer_arrowhead_offset(
                m_start_item, end_qcoords, m_start_item.getAtom(),
                m_end_item.getAtom(), m_is_secondary_connection);
        }
        if (end_has_arrowhead) {
            end_offset = get_monomer_arrowhead_offset(
                m_end_item, start_qcoords, m_end_item.getAtom(),
                m_start_item.getAtom(), m_is_secondary_connection);
        }
    } else {
        const auto start_magnitude = m_start_item.boundingRect().height() / 2 +
                                     MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
        const auto end_magnitude = m_end_item.boundingRect().height() / 2 +
                                   MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
        // Positive lanes are above the backbone. Qt's y axis points down, so
        // an above-backbone offset has a negative y component.
        const auto vertical_direction = m_connector_lane > 0 ? -1 : 1;
        start_offset.ry() = vertical_direction * start_magnitude;
        end_offset.ry() = vertical_direction * end_magnitude;
    }

    auto end_pos = end_qcoords - start_qcoords + end_offset;

    // Every numbered lane floats away from its endpoint anchors. This keeps a
    // horizontal line from running through a diamond on another connection.
    // Larger lane numbers move outward by one additional diamond width.
    qreal lane_offset_y = 0;
    if (m_connector_lane != 0) {
        const int lane_distance = std::abs(m_connector_lane);
        const auto vertical_direction = m_connector_lane > 0 ? -1 : 1;
        lane_offset_y = vertical_direction * lane_distance * 2 *
                        MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
    }
    const QPointF line_displacement(0, lane_offset_y);
    m_connector_line =
        QLineF(start_offset + line_displacement, end_pos + line_displacement);

    if (m_connector_lane != 0) {
        // Complete the routed connector with one line from each residue to the
        // horizontal lane. Together, the three lines form two right-angle
        // bends. Diamond connections begin at the diamond's outward tip;
        // connections without diamonds begin at the residue edge.
        const auto vertical_direction = m_connector_lane > 0 ? -1 : 1;
        const auto get_join_start =
            [vertical_direction](const QPointF& endpoint,
                                 const bool has_diamond) {
                const auto offset_direction =
                    has_diamond ? vertical_direction : -vertical_direction;
                return endpoint +
                       QPointF(0, offset_direction *
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
            };
        m_start_join_line =
            QLineF(get_join_start(start_offset, start_has_arrowhead),
                   m_connector_line.p1());
        m_end_join_line = QLineF(get_join_start(end_pos, end_has_arrowhead),
                                 m_connector_line.p2());
    }
    if (start_has_arrowhead) {
        add_diamond_arrowhead_to_path(m_arrowhead_path, start_offset,
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                      offset_is_at_corner(start_offset));
    }
    if (end_has_arrowhead) {
        add_diamond_arrowhead_to_path(m_arrowhead_path, end_pos,
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                      offset_is_at_corner(end_offset));
    }
    m_midpoint = m_connector_line.center();
    m_selection_highlighting_path = path_around_line(
        m_connector_line, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH);
    m_predictive_highlighting_path = path_around_line(
        m_connector_line, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH);
    for (const auto& join_line : {m_start_join_line, m_end_join_line}) {
        if (!join_line.isNull()) {
            m_selection_highlighting_path |= path_around_line(
                join_line, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH);
            m_predictive_highlighting_path |= path_around_line(
                join_line, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH);
        }
    }
    if (start_has_arrowhead) {
        or_diamond_arrowhead_to_path(m_selection_highlighting_path,
                                     start_offset,
                                     BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                     offset_is_at_corner(start_offset));
        or_diamond_arrowhead_to_path(m_predictive_highlighting_path,
                                     start_offset,
                                     BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                     offset_is_at_corner(start_offset));
    }
    if (end_has_arrowhead) {
        or_diamond_arrowhead_to_path(m_selection_highlighting_path, end_pos,
                                     BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                     offset_is_at_corner(end_offset));
        or_diamond_arrowhead_to_path(m_predictive_highlighting_path, end_pos,
                                     BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH +
                                         MONOMER_CONNECTOR_ARROWHEAD_RADIUS,
                                     offset_is_at_corner(end_offset));
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
    if (!m_start_join_line.isNull()) {
        painter->drawLine(m_start_join_line);
    }
    if (!m_end_join_line.isNull()) {
        painter->drawLine(m_end_join_line);
    }

    if (!m_arrowhead_path.isEmpty()) {
        painter->setPen(m_arrowhead_pen);
        painter->setBrush(m_arrowhead_brush);
        painter->drawPath(m_arrowhead_path);
    }

    painter->restore();
}

QColor MonomerConnectorItem::getConnectorColor() const
{
    return m_is_dark_mode ? m_connector_color_dark_bg : m_connector_color;
}

void MonomerConnectorItem::setConnectorStyle(const QColor& connector_color,
                                             const qreal connector_width)
{
    m_connector_pen.setColor(connector_color);
    m_connector_pen.setWidthF(connector_width);
    m_arrowhead_pen.setColor(connector_color);
    m_arrowhead_pen.setWidthF(connector_width);
    m_arrowhead_brush.setColor(connector_color);
}

} // namespace sketcher
} // namespace schrodinger
