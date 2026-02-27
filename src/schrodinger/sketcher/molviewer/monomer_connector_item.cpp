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
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

std::unordered_map<ConnectorType, std::tuple<QColor, qreal, Qt::PenStyle>>
    PEN_STYLE_FOR_CONNECTOR_TYPE = {
        {ConnectorType::CHEM,
         {CHEM_CONNECTOR_COLOR, CHEM_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_LINEAR,
         {AA_LINEAR_CONNECTOR_COLOR, AA_LINEAR_CONNECTOR_WIDTH,
          Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_BRANCHING,
         {AA_BRANCHING_CONNECTOR_COLOR, AA_BRANCHING_CONNECTOR_WIDTH,
          Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_DISULFIDE,
         {DISULFIDE_CONNECTOR_COLOR, DISULFIDE_CONNECTOR_WIDTH,
          Qt::PenStyle::SolidLine}},
        {ConnectorType::PEPTIDE_SIDE_CHAIN,
         {AA_LINEAR_CONNECTOR_COLOR, AA_LINEAR_CONNECTOR_WIDTH,
          Qt::PenStyle::SolidLine}},
        {ConnectorType::NA_BASE,
         {NA_BASE_CONNECTOR_COLOR, NA_BASE_CONNECTOR_WIDTH,
          Qt::PenStyle::DotLine}},
        {ConnectorType::NA_BACKBONE,
         {NA_BACKBONE_CONNECTOR_COLOR, NA_BACKBONE_CONNECTOR_WIDTH,
          Qt::PenStyle::SolidLine}},
        {ConnectorType::NA_BACKBONE_TO_BASE,
         {NA_BACKBONE_TO_BASE_CONNECTOR_COLOR,
          NA_BACKBONE_TO_BASE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine}}};

// clang-format off
const std::unordered_map<ConnectorType, QColor> DARK_BG_COLOR_FOR_CONNECTOR_TYPE = {
    {ConnectorType::CHEM,                CHEM_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::PEPTIDE_LINEAR,      AA_LINEAR_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::PEPTIDE_BRANCHING,   AA_BRANCHING_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::PEPTIDE_DISULFIDE,   DISULFIDE_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::PEPTIDE_SIDE_CHAIN,  AA_LINEAR_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::NA_BASE,             NA_BASE_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::NA_BACKBONE,         NA_BACKBONE_CONNECTOR_COLOR_DARK_BG},
    {ConnectorType::NA_BACKBONE_TO_BASE, NA_BACKBONE_TO_BASE_CONNECTOR_COLOR_DARK_BG}};
// clang-format on

} // namespace

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
    auto connector_type = get_connector_type(m_bond, m_is_secondary_connection);
    auto [start_has_arrowhead, end_has_arrowhead] =
        does_connector_have_arrowheads(m_bond, connector_type);
    auto [connector_color, connector_width, connector_pen_style] =
        PEN_STYLE_FOR_CONNECTOR_TYPE.at(connector_type);
    m_connector_color = connector_color;
    m_connector_color_dark_bg =
        DARK_BG_COLOR_FOR_CONNECTOR_TYPE.at(connector_type);
    auto color = getConnectorColor();
    m_connector_pen = QPen(color, connector_width, connector_pen_style);
    m_arrowhead_pen = QPen(color, connector_width);
    m_arrowhead_pen.setJoinStyle(Qt::PenJoinStyle::MiterJoin);
    m_arrowhead_brush = QBrush(color);

    auto start_qcoords = m_start_item.pos();
    auto end_qcoords = m_end_item.pos();
    setPos(start_qcoords);

    QPointF start_offset;
    if (start_has_arrowhead) {
        start_offset.ry() -=
            get_monomer_arrowhead_offset(m_start_item, end_qcoords);
        add_diamond_arrowhead_to_path(m_arrowhead_path, start_offset,
                                      MONOMER_CONNECTOR_ARROWHEAD_RADIUS);
    }

    auto end_pos = end_qcoords - start_qcoords;
    if (end_has_arrowhead) {
        end_pos.ry() -= get_monomer_arrowhead_offset(m_end_item, start_qcoords);
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

void MonomerConnectorItem::setDarkMode(bool is_dark)
{
    m_is_dark_mode = is_dark;
    auto color = getConnectorColor();
    m_connector_pen.setColor(color);
    m_arrowhead_pen.setColor(color);
    m_arrowhead_brush.setColor(color);
    update();
}

QColor MonomerConnectorItem::getConnectorColor() const
{
    return m_is_dark_mode ? m_connector_color_dark_bg : m_connector_color;
}

} // namespace sketcher
} // namespace schrodinger
