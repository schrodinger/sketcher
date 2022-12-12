#include "schrodinger/sketcher/molviewer/bond_item.h"

#include <GraphMol/ROMol.h>

#include <QtGlobal>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>

#include "schrodinger/sketcher/molviewer/atom_item.h"

namespace schrodinger
{
namespace sketcher
{

BondItem::BondItem(RDKit::Bond* bond, const AtomItem& start_item,
                   const AtomItem& end_item, BondItemSettings& settings,
                   QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_bond(bond),
    m_start_item(start_item),
    m_end_item(end_item),
    m_settings(settings)
{
    m_solid_pen.setJoinStyle(Qt::RoundJoin);
    m_solid_pen.setCapStyle(Qt::RoundCap);
    m_dashed_pen.setJoinStyle(Qt::RoundJoin);
    m_dashed_pen.setCapStyle(Qt::RoundCap);
    m_dashed_pen.setDashPattern({3.0, 3.0});

    updateCachedData();
}

int BondItem::type() const
{
    return Type;
}

void BondItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();

    m_solid_pen.setWidthF(m_settings.m_bond_width);
    m_dashed_pen.setWidthF(m_settings.m_bond_width);

    // move this item so it's on top of the start atom
    setPos(m_start_item.pos());
    // the bond line and bounding rect must both be relative to the start atom
    // since we just set that as the pos of this item
    QPointF bond_end = m_end_item.pos() - m_start_item.pos();
    QLineF bond_line = QLineF(QPointF(0, 0), bond_end);
    m_to_paint = calculateLinesToPaint(bond_line);
    m_shape = pathAroundLine(bond_line, PREDICTIVE_HIGHLIGHTING_HALF_WIDTH);
    m_bounding_rect = m_shape.boundingRect();

    // TODO: calculate bond intersections, store clipping path
    // TODO: deal with atom radii for aromatics
}

std::vector<ToPaint>
BondItem::calculateLinesToPaint(const QLineF& bond_line) const
{
    RDKit::Bond::BondType bond_type = m_bond->getBondType();
    RDKit::Bond::BondDir bond_direction = m_bond->getBondDir();

    QLineF trimmed_line(bond_line);
    trimLineToBoundAtoms(trimmed_line);

    if (bond_type == RDKit::Bond::BondType::SINGLE &&
        bond_direction == RDKit::Bond::BondDir::NONE) {
        return {{m_solid_pen, {trimmed_line}, {}}};
    } else if (bond_type == RDKit::Bond::BondType::ZERO) {
        return {{m_dashed_pen, {trimmed_line}, {}}};
    } else if (bond_type == RDKit::Bond::BondType::TRIPLE) {
        auto bond_lines = calcTripleBondLines(trimmed_line);
        return {{m_solid_pen, bond_lines, {}}};
    } else {
        return {};
    }
}

std::vector<QLineF>
BondItem::calcTripleBondLines(const QLineF& trimmed_line) const
{
    QLineF normal = trimmed_line.normalVector();
    normal.setLength(m_settings.m_double_bond_spacing);
    QPointF offset = normal.p2() - normal.p1();
    QLineF offset_line1 = trimmed_line.translated(offset);
    QLineF offset_line2 = trimmed_line.translated(-offset);
    return {offset_line1, trimmed_line, offset_line2};
}

void BondItem::trimLineToBoundAtoms(QLineF& line) const
{
    // trim to the start atom
    for (const QRectF& subrect : m_start_item.getSubrects()) {
        // this bond uses the same local coordinate system as the start
        // atom, so we don't need to map these subrects
        trimLineToRect(line, subrect);
    }
    // trim to the end atom
    for (const QRectF& subrect : m_end_item.getSubrects()) {
        QRectF mapped_rect = mapRectFromItem(&m_end_item, subrect);
        trimLineToRect(line, mapped_rect);
    }
}

void BondItem::trimLineToRect(QLineF& line, const QRectF& rect) const
{
    // expand rect by 4 pixels in all directions
    const qreal margin_size = 4.0;
    const QMarginsF margins(margin_size, margin_size, margin_size, margin_size);
    const QRectF enlarged_rect = rect + margins;

    QPointF inter_point; // the intersection point
    if (enlarged_rect.contains(line.p1()) &&
        intersectionOfLineAndRect(line, enlarged_rect, inter_point)) {
        line.setP1(inter_point);
    } else if (enlarged_rect.contains(line.p2()) &&
               intersectionOfLineAndRect(line, enlarged_rect, inter_point)) {
        line.setP2(inter_point);
    }
}

bool BondItem::intersectionOfLineAndRect(const QLineF& line, const QRectF& rect,
                                         QPointF& inter_point) const
{
    // break the bounding rect up into four lines
    QLineF top(rect.topLeft(), rect.topRight());
    QLineF left(rect.topLeft(), rect.bottomLeft());
    QLineF bottom(rect.bottomLeft(), rect.bottomRight());
    QLineF right(rect.topRight(), rect.bottomRight());

    bool found_intersection =
        (line.intersects(top, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(left, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(bottom, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(right, &inter_point) == QLineF::BoundedIntersection);
    return found_intersection;
}

QPainterPath BondItem::pathAroundLine(const QLineF& line,
                                      const qreal half_width) const
{
    QLineF normal = line.normalVector();
    normal.setLength(half_width);
    QPointF offset = normal.p2() - normal.p1();
    QPointF p1 = line.p1();
    QPointF p2 = line.p2();
    QPainterPath path;
    path.moveTo(p1 + offset);
    path.lineTo(p2 + offset);
    path.lineTo(p2 - offset);
    path.lineTo(p1 - offset);
    path.closeSubpath();
    return path;
}

QPainterPath BondItem::shape() const
{
    return QPainterPath(m_shape);
}

void BondItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget)
{
    painter->save();
    for (auto to_paint : m_to_paint) {
        painter->setPen(to_paint.pen);
        for (auto line : to_paint.lines) {
            painter->drawLine(line);
        }
        for (auto polygon : to_paint.polygons) {
            painter->drawPolygon(polygon);
        }
    }
    painter->restore();
}

} // namespace sketcher
} // namespace schrodinger
