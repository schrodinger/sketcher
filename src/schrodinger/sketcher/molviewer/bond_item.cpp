#include "schrodinger/sketcher/molviewer/bond_item.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_set>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/QueryOps.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RingInfo.h>
#include <QMarginsF>
#include <QPainter>
#include <QPointF>
#include <QTransform>
#include <Qt>
#include <QtGlobal>
#include <QtMath>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"

namespace schrodinger
{
namespace sketcher
{

struct NumBondsInRing {
    unsigned int num_bonds;

    // Double bonds and aromatic bonds are drawn with an asymmetric second line
    // that needs to be placed inside the ring. This is not true for any other
    // bond, including the double/aromatic bond query

    unsigned int num_double_or_aromatic_bonds;
    int ring_index;

    // See findBestRingForBond documentation for an explanation of the
    // ordering criteria
    bool operator<(const NumBondsInRing& other) const
    {
        const unsigned int& cutoff = DOUBLE_BOND_BEST_RING_SIZE_CUTOFF;
        if (num_bonds <= cutoff && other.num_bonds > cutoff) {
            return false;
        }
        if (num_bonds > cutoff && other.num_bonds <= cutoff) {
            return true;
        }
        if (num_double_or_aromatic_bonds !=
            other.num_double_or_aromatic_bonds) {
            return num_double_or_aromatic_bonds <
                   other.num_double_or_aromatic_bonds;
        }
        return ring_index < other.ring_index;
    };
};

BondItem::BondItem(const RDKit::Bond* bond, const AtomItem& start_item,
                   const AtomItem& end_item, const Fonts& fonts,
                   BondItemSettings& settings, QGraphicsItem* parent) :

    AbstractGraphicsItem(parent),
    m_bond(bond),
    m_start_item(start_item),
    m_end_item(end_item),
    m_fonts(fonts),
    m_settings(settings)
{
    m_solid_pen.setJoinStyle(Qt::RoundJoin);
    m_solid_pen.setCapStyle(Qt::RoundCap);
    m_dashed_pen.setJoinStyle(Qt::RoundJoin);
    m_dashed_pen.setCapStyle(Qt::RoundCap);
    m_dashed_pen.setDashPattern({3.0, 3.0});
    m_chirality_pen = QPen(CHIRALITY_LABEL_COLOR);

    setZValue(static_cast<qreal>(ZOrder::BOND));
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

    m_solid_pen.setColor(m_settings.m_color);
    m_dashed_pen.setColor(m_settings.m_color);
    m_solid_brush.setColor(m_settings.m_color);

    m_solid_pen.setWidthF(m_settings.m_bond_width);
    m_dashed_pen.setWidthF(m_settings.m_bond_width);

    // move this item so it's on top of the start atom
    setPos(m_start_item.pos());
    // the bond line and bounding rect must both be relative to the start atom
    // since we just set that as the pos of this item
    QPointF bond_end = m_end_item.pos() - m_start_item.pos();
    QLineF bond_line = QLineF(QPointF(0, 0), bond_end);
    m_to_paint = calculateLinesToPaint(bond_line);
    m_selection_highlighting_path =
        pathAroundLine(bond_line, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH);
    m_predictive_highlighting_path =
        pathAroundLine(bond_line, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH);
    m_shape = QPainterPath(m_predictive_highlighting_path);
    m_bounding_rect = m_shape.boundingRect();

    QString label_text;
    bool is_query = m_bond->hasQuery();
    if (is_query) {
        label_text = QString::fromStdString(m_bond->getQuery()->getTypeLabel());
    } else {
        label_text = get_bond_stereo_label(*m_bond);
    }

    // add the label (if present) to the bounding rect
    if (!label_text.isEmpty()) {

        // TODO: make bond partially transparent behind query label
        std::tie(m_text_angle, m_text_pos, m_text_size) =
            getStereoAnnotationParameters(label_text, !is_query);
        auto rect = QRectF(m_text_pos, m_text_size);
        rect.moveCenter(QPointF(0, 0));
        auto rotation = QTransform().rotate(m_text_angle);
        auto bounding_rect = rotation.mapRect(rect);

        bounding_rect.moveCenter(m_text_pos);
        m_bounding_rect = m_bounding_rect.united(bounding_rect);
    }
    m_label_text = label_text;

    // TODO: calculate bond intersections, store clipping path
    // TODO: deal with atom radii for aromatics
}

std::vector<ToPaint>
BondItem::calculateLinesToPaint(const QLineF& bond_line) const
{
    RDKit::Bond::BondType bond_type = m_bond->getBondType();
    RDKit::Bond::BondDir bond_direction = m_bond->getBondDir();
    QLineF trimmed_line = trimLineToBoundAtoms(bond_line);

    switch (bond_type) {
        case RDKit::Bond::BondType::SINGLE:
            switch (bond_direction) {
                case RDKit::Bond::BondDir::BEGINWEDGE: {
                    // an "up" chiral bond is drawn with a wedge
                    QPolygonF wedge = calcSolidWedge(trimmed_line);
                    return {{m_solid_pen, {}, {wedge}}};
                }
                case RDKit::Bond::BondDir::BEGINDASH: {
                    // a "down" chiral bond is drawn with a dashed wedge
                    auto dashes = calcDashedWedge(trimmed_line);
                    return {{m_solid_pen, dashes, {}}};
                }
                case RDKit::Bond::BondDir::UNKNOWN: {
                    // an "up or down" bond is drawn with a zig-zag wedge line
                    auto lines = calcSquigglyWedge(trimmed_line);
                    return {{m_solid_pen, lines, {}}};
                }
                default:
                    // normal non-chiral single bond
                    return {{m_solid_pen, {trimmed_line}, {}}};
            }
        case RDKit::Bond::BondType::DATIVE:
        case RDKit::Bond::BondType::DATIVEONE:
        case RDKit::Bond::BondType::DATIVEL:
        case RDKit::Bond::BondType::DATIVER: {
            // a dative bond, a.k.a. a coordinate bond, is drawn with an arrow
            QPolygonF arrow_tip = calcArrowTip(trimmed_line);
            return {{m_solid_pen, {trimmed_line}, {arrow_tip}}};
        }
        case RDKit::Bond::BondType::ZERO:
            return {{m_dashed_pen, {trimmed_line}, {}}};
        case RDKit::Bond::BondType::DOUBLE: {
            QLineF line_one, line_two;
            if (bond_direction == RDKit::Bond::BondDir::EITHERDOUBLE) {
                // a "cis or trans" double bond is drawn with two crossed lines
                std::tie(line_one, line_two) =
                    calcSymmetricDoubleBondLines(trimmed_line);
                crossDoubleBondLines(line_one, line_two);
            } else {
                std::tie(line_one, line_two) =
                    calcDoubleBondLines(trimmed_line, bond_line);
            }
            return {{m_solid_pen, {line_one, line_two}, {}}};
        }
        case RDKit::Bond::BondType::AROMATIC: {
            auto [line_one, line_two] =
                calcDoubleBondLines(trimmed_line, bond_line);
            return {{m_solid_pen, {line_one}, {}},
                    {m_dashed_pen, {line_two}, {}}};
        }
        case RDKit::Bond::BondType::TRIPLE: {
            auto bond_lines = calcTripleBondLines(trimmed_line);
            return {{m_solid_pen, bond_lines, {}}};
        }
        default:
            // unrecognized bond type, so paint it as a single bond
            // TODO: report something here?  See SKETCH-1852.
            return {{m_solid_pen, {trimmed_line}, {}}};
    }
}

QPolygonF BondItem::calcArrowTip(const QLineF& trimmed_line) const
{
    QLineF parallel_vector = trimmed_line.unitVector();
    QLineF normal_vector = parallel_vector.normalVector();
    parallel_vector.setLength(DATIVE_ARROW_LENGTH);
    QPointF arrow_length_offset = parallel_vector.p1() - parallel_vector.p2();
    normal_vector.setLength(DATIVE_ARROW_HALF_WIDTH);
    QPointF arrow_half_width_offset = normal_vector.p2() - normal_vector.p1();
    QPointF line_end = trimmed_line.p2();

    QPointF arrow_point_1 =
        line_end + arrow_length_offset + arrow_half_width_offset;
    QPointF arrow_point_2 =
        line_end + arrow_length_offset - arrow_half_width_offset;
    QPolygonF arrow_tip({line_end, arrow_point_1, arrow_point_2, line_end});
    return arrow_tip;
}

std::tuple<QPointF, QPointF, QPointF>
BondItem::calcWedgeEnds(const QLineF& trimmed_line) const
{
    QPointF wedge_start = trimmed_line.p1();
    QPointF line_end = trimmed_line.p2();
    QLineF normal_vector = trimmed_line.normalVector();
    normal_vector.setLength(m_settings.m_wedge_width / 2.0);
    QPointF offset = normal_vector.p2() - normal_vector.p1();
    return std::make_tuple(wedge_start, line_end + offset, line_end - offset);
}

QPolygonF BondItem::calcSolidWedge(const QLineF& trimmed_line) const
{
    auto [wedge_start, p2, p3] = calcWedgeEnds(trimmed_line);
    return QPolygonF({wedge_start, p2, p3, wedge_start});
}

std::tuple<unsigned int, QPointF, QPointF, QPointF, std::vector<QLineF>>
BondItem::calcDashedWedgeParams(const QLineF& trimmed_line) const
{
    auto [wedge_start, p2, p3] = calcWedgeEnds(trimmed_line);
    unsigned int num_dashes = trimmed_line.length() / m_settings.m_hash_spacing;
    num_dashes = std::clamp(num_dashes, 0u, 100u);
    QPointF offset_towards_p2 = (p2 - wedge_start) / num_dashes;
    QPointF offset_towards_p3 = (p3 - wedge_start) / num_dashes;
    std::vector<QLineF> dashes;
    dashes.reserve(num_dashes + 1);
    return std::make_tuple(num_dashes, wedge_start, offset_towards_p2,
                           offset_towards_p3, dashes);
}

std::vector<QLineF> BondItem::calcDashedWedge(const QLineF& trimmed_line) const
{
    auto [num_dashes, wedge_start, offset_towards_p2, offset_towards_p3,
          dashes] = calcDashedWedgeParams(trimmed_line);
    // Note that we intentionally use <= here and < in calcSquigglyWedge.  In
    // this method, we want one last dash at the end of trimmed_line.  In
    // calcSquigglyWedge, we don't want a squiggle starting at the end of
    // trimmed_line and then going m_hash_spacing pixels closer to the atom.
    for (unsigned int i = 0; i <= num_dashes; ++i) {
        // draw lines across the wedge
        dashes.emplace_back(wedge_start + offset_towards_p2 * i,
                            wedge_start + offset_towards_p3 * i);
    }
    return dashes;
}

std::vector<QLineF>
BondItem::calcSquigglyWedge(const QLineF& trimmed_line) const
{
    auto [num_dashes, wedge_start, offset_towards_p2, offset_towards_p3,
          dashes] = calcDashedWedgeParams(trimmed_line);
    for (unsigned int i = 0; i < num_dashes; ++i) {
        // alternate between zigging and zagging across the wedge
        unsigned int is_odd = i % 2;
        unsigned int is_even = 1 - is_odd;
        dashes.emplace_back(wedge_start + offset_towards_p2 * (i + is_even),
                            wedge_start + offset_towards_p3 * (i + is_odd));
    }
    return dashes;
}

std::pair<QLineF, QLineF>
BondItem::calcDoubleBondLines(const QLineF& trimmed_line,
                              const QLineF& bond_line) const
{
    if (isDoubleBondSymmetric()) {
        return calcSymmetricDoubleBondLines(trimmed_line);
    } else {
        return calcAsymmetricDoubleBondLines(trimmed_line, bond_line);
    }
}

bool BondItem::isDoubleBondSymmetric() const
{
    RDKit::Atom* begin_atom = m_bond->getBeginAtom();
    RDKit::Atom* end_atom = m_bond->getEndAtom();
    int begin_degree = begin_atom->getDegree();
    int end_degree = end_atom->getDegree();
    if (begin_degree == 1 && end_degree == 1) {
        // We always draw the bond as symmetrical when neither atom has any
        // other bonds (i.e. this bond is completely isolated)
        return true;
    }
    if ((begin_degree == 1 &&
         (end_degree == 3 || RDKit::queryIsAtomInRing(end_atom))) ||
        (end_degree == 1 &&
         (begin_degree == 3 || RDKit::queryIsAtomInRing(begin_atom)))) {
        // We always draw the bond as symmetrical when one atom has no other
        // bonds (i.e. the bond is terminal) and the other atom has exactly two
        // other bonds (or it's part of a ring, in which case it's allowed to
        // have more than two)
        return true;
    }
    if ((begin_degree == 1 || m_start_item.labelIsVisible()) &&
        (end_degree == 1 || m_end_item.labelIsVisible()) &&
        !RDKit::queryIsBondInRing(m_bond)) {
        // We always draw the bond as symmetrical when this bond is not part of
        // a ring and the labels for both atoms are either visible or terminal
        // atoms
        return true;
    }
    // "cis or trans" bonds (i.e. crossed double bonds) are always drawn
    // symmetrically, but calculateLinesToPaint calls
    // calcSymmetricDoubleBondLines directly in that scenario, so we don't have
    // to worry about checking that here.

    // The bond should be drawn asymmetrically
    return false;
}

std::pair<QLineF, QLineF>
BondItem::calcSymmetricDoubleBondLines(const QLineF& trimmed_line) const
{
    qreal offset = m_settings.m_double_bond_spacing / 2.0;
    auto [offset_line_one, offset_line_two] =
        calcOffsetBondLines(trimmed_line, offset);
    return std::make_pair(offset_line_one, offset_line_two);
}

std::pair<QLineF, QLineF>
BondItem::calcAsymmetricDoubleBondLines(const QLineF& trimmed_line,
                                        const QLineF& bond_line) const
{
    QPointF offset = calcAsymmetricDoubleBondOffset(trimmed_line);
    QLineF inner_line = trimmed_line.translated(offset);
    trimDoubleBondInnerLine(inner_line, offset, bond_line);
    return std::make_pair(trimmed_line, inner_line);
}

QPointF
BondItem::calcAsymmetricDoubleBondOffset(const QLineF& trimmed_line) const
{
    // Create an offset of the appropriate length
    QLineF normal = trimmed_line.normalVector();
    normal.setLength(m_settings.m_double_bond_spacing);
    QPointF offset = normal.p2() - normal.p1();

    // Figure out whether offset is on the correct side of the line, i.e.,
    // should we be drawing the second line to the left or the right of
    // trimmed_line.  If offset is on the wrong side, flip it.
    const RDKit::ROMol& molecule = m_bond->getOwningMol();
    RDKit::RingInfo* ring_info = molecule.getRingInfo();
    int ring_index = findBestRingForBond(molecule, ring_info);
    if (ring_index >= 0) {
        // this bond is in a ring, so we want to draw the second line *inside*
        // of the ring
        QPointF ring_center = findRingCenter(molecule, ring_info, ring_index);
        // make sure offset is on the same side of the bond as ring_center
        if (!are_points_on_same_side_of_line(offset, ring_center,
                                             trimmed_line.p2())) {
            // offset is on the wrong side of the bond, so flip it
            offset *= -1;
        }
    } else {
        // this bond isn't in a ring, so we want the second line to be on the
        // side of the bond with more neighboring atoms.
        int num_same_side = 0;
        int num_opposite_side = 0;
        const RDKit::Atom* begin_atom = m_bond->getBeginAtom();
        const RDKit::Atom* end_atom = m_bond->getEndAtom();
        const RDKit::Conformer& conf = molecule.getConformer();
        QPointF line_endpoint = trimmed_line.p2();
        // count the atoms that neighbor begin_atom (other than end_atom)
        countNeighboringAtomsBySide(molecule, begin_atom, end_atom, conf,
                                    offset, line_endpoint, num_same_side,
                                    num_opposite_side);
        // count the atoms that neighbor end_atom (other than begin_atom)
        countNeighboringAtomsBySide(molecule, end_atom, begin_atom, conf,
                                    offset, line_endpoint, num_same_side,
                                    num_opposite_side);
        if (num_opposite_side > num_same_side) {
            // offset is on the wrong side of the bond, so flip it
            offset *= -1;
        } else if (num_opposite_side == num_same_side && offset.y() > 0) {
            // if there's a tie, always put the second line below the first one
            offset *= -1;
        }
    }
    return offset;
}

int BondItem::findBestRingForBond(const RDKit::ROMol& molecule,
                                  const RDKit::RingInfo* ring_info) const
{
    // figure out what rings this bond is a part of
    std::vector<int> ring_indices = ring_info->bondMembers(m_bond->getIdx());
    if (ring_indices.empty()) {
        // this bond isn't part of any rings
        return -1;
    }
    if (ring_indices.size() == 1) {
        // this bond is part of only one ring
        return ring_indices[0];
    }

    // This bond is part of multiple rings, so we need to pick the "best" one.
    // First, we count the number of bonds and double bonds in each ring, since
    // that's what we use to decide.
    std::vector<NumBondsInRing> num_bonds_in_rings;
    num_bonds_in_rings.reserve(ring_indices.size());
    auto bond_rings = ring_info->bondRings();
    for (int ring_index : ring_indices) {
        auto bond_indices = bond_rings[ring_index];
        unsigned int num_bonds = bond_indices.size();
        unsigned int num_double_or_aromatic_bonds = 0;
        for (int bond_index : bond_indices) {
            const RDKit::Bond* cur_bond = molecule.getBondWithIdx(bond_index);
            if (cur_bond->getBondType() == RDKit::Bond::BondType::DOUBLE ||
                cur_bond->getBondType() == RDKit::Bond::BondType::AROMATIC) {
                ++num_double_or_aromatic_bonds;
            }
        }
        num_bonds_in_rings.push_back(
            {num_bonds, num_double_or_aromatic_bonds, ring_index});
    }

    // Now pick the "best" ring, i.e., the ring we should draw this bond inside
    // of.  See the header file for documentation of the selection criteria.
    NumBondsInRing best_ring =
        *std::max_element(num_bonds_in_rings.begin(), num_bonds_in_rings.end());
    return best_ring.ring_index;
}

QPointF BondItem::findRingCenter(const RDKit::ROMol& molecule,
                                 const RDKit::RingInfo* ring_info,
                                 const int ring_index) const
{
    auto atom_rings = ring_info->atomRings();
    const auto& conformer = molecule.getConformer();
    RDGeom::Point3D center;
    auto atom_indices = atom_rings[ring_index];
    for (int atom_index : atom_indices) {
        auto pos = conformer.getAtomPos(atom_index);
        center += pos;
    }
    center /= atom_indices.size();
    return mapFromScene(to_scene_xy(center));
}

void BondItem::countNeighboringAtomsBySide(
    const RDKit::ROMol& molecule, const RDKit::Atom* neighbors_of,
    const RDKit::Atom* atom_to_skip, const RDKit::Conformer& conf,
    const QPointF& offset, const QPointF& line_endpoint, int& num_same_side,
    int& num_opposite_side) const
{
    for (const RDKit::Atom* neighbor : molecule.atomNeighbors(neighbors_of)) {
        if (neighbor == atom_to_skip) {
            continue;
        }
        RDGeom::Point3D neighbor_coords = conf.getAtomPos(neighbor->getIdx());
        QPointF neighbor_qcoords = mapFromScene(to_scene_xy(neighbor_coords));
        if (are_points_on_same_side_of_line(neighbor_qcoords, offset,
                                            line_endpoint)) {
            ++num_same_side;
        } else {
            ++num_opposite_side;
        }
    }
}

void BondItem::trimDoubleBondInnerLine(QLineF& inner_line,
                                       const QPointF& offset,
                                       const QLineF& bond_line) const
{
    QLineF shortening_line = QLineF(inner_line);
    shortening_line.setLength(DOUBLE_BOND_INNER_LINE_SHORTENING);
    QPointF shortening = shortening_line.p2() - shortening_line.p1();

    int begin_degree = m_bond->getBeginAtom()->getDegree();
    // note that we're working in local coordinates here, so bond_line.p1() is
    // (0, 0)
    if (begin_degree > 1 && QLineF(offset, inner_line.p1()).length() <
                                DOUBLE_BOND_INNER_LINE_SHORTENING) {
        // the starting atom side of the bond hasn't been trimmed enough, so
        // trim it by 6 pixels
        inner_line.setP1(offset + shortening);
    }

    int end_degree = m_bond->getEndAtom()->getDegree();
    QPointF end_offset = bond_line.p2() + offset;
    if (end_degree > 1 && QLineF(end_offset, inner_line.p2()).length() <
                              DOUBLE_BOND_INNER_LINE_SHORTENING) {
        // the ending atom side of the bond hasn't been trimmed enough, so trim
        // it by 6 pixels
        inner_line.setP2(end_offset - shortening);
    }
}

void BondItem::crossDoubleBondLines(QLineF& line_one, QLineF& line_two) const
{
    QPointF end_one = line_one.p2();
    QPointF end_two = line_two.p2();
    line_one.setP2(end_two);
    line_two.setP2(end_one);
}

std::pair<QLineF, QLineF>
BondItem::calcOffsetBondLines(const QLineF& center_line,
                              const qreal distance) const
{
    QLineF normal = center_line.normalVector();
    normal.setLength(distance);
    QPointF offset = normal.p2() - normal.p1();
    return std::make_pair(center_line.translated(offset),
                          center_line.translated(-offset));
}

std::vector<QLineF>
BondItem::calcTripleBondLines(const QLineF& trimmed_line) const
{
    auto [offset_line_one, offset_line_two] =
        calcOffsetBondLines(trimmed_line, m_settings.m_double_bond_spacing);
    return {offset_line_one, trimmed_line, offset_line_two};
}

QLineF BondItem::trimLineToBoundAtoms(const QLineF& line) const
{
    QLineF trimmed_line(line);
    // trim to the start atom
    for (const QRectF& subrect : m_start_item.getSubrects()) {
        // this bond uses the same local coordinate system as the start
        // atom, so we don't need to map these subrects
        trim_line_to_rect(trimmed_line, subrect);
    }
    // trim to the end atom
    for (const QRectF& subrect : m_end_item.getSubrects()) {
        QRectF mapped_rect = mapRectFromItem(&m_end_item, subrect);
        trim_line_to_rect(trimmed_line, mapped_rect);
    }
    return trimmed_line;
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

void BondItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget)
{
    painter->save();
    painter->setBrush(m_solid_brush);
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
    if (!m_label_text.isEmpty()) {
        painter->save();
        painter->setPen(m_chirality_pen);
        painter->setFont(m_fonts.m_chirality_font);
        paintAnnotation(painter, m_text_angle, m_text_pos, m_text_size,
                        m_label_text);
        painter->restore();
    }
}

const RDKit::Bond* BondItem::getBond() const
{
    return m_bond;
}

std::tuple<qreal, QPointF, QSizeF>
BondItem::getStereoAnnotationParameters(const QString& label,
                                        const bool draw_text_above_bond) const
{
    auto distance = draw_text_above_bond
                        ? VIEW_SCALE * BOND_STEREO_LABEL_DISTANCE_RATIO
                        : 0;
    auto bond_end = m_end_item.pos() - m_start_item.pos();
    auto bond_line = QLineF(QPointF(0, 0), bond_end);
    // Calculate the unit vector in the direction of the line segment
    auto unit_vector_line = bond_line.unitVector();

    // Calculate the offset vector perpendicular to the unit vector, always
    // pointing towards the upper hemiplane (which is y < 0, since
    // QGraphicsScenes's y axis grows downwards)
    auto normal_line = unit_vector_line.normalVector();
    auto offset_vector = normal_line.p2() * distance;
    if (offset_vector.y() > 0) {
        offset_vector = -offset_vector;
    }

    // Calculate the position of the text as the midpoint between the two points
    // plus the offset
    auto text_pos = bond_end * 0.5 + offset_vector;

    // Calculate the angle between the x-axis and the line segment in degrees so
    // that the angle is between -90 and 90
    auto angle = std::remainder(bond_line.angle(), 180.0);
    QFont font = m_fonts.m_chirality_font;
    auto text_size = QFontMetrics(font).size(Qt::TextSingleLine, label);
    return std::make_tuple(angle, text_pos, text_size);
}

void BondItem::paintAnnotation(QPainter* painter, qreal angle,
                               const QPointF& text_pos, const QSizeF& text_size,
                               const QString& text)
{
    painter->save();
    // Position and rotate the painter to the correct angle
    painter->translate(text_pos);
    painter->rotate(-angle);

    // Draw the text centered at the rotated position
    QFont font = painter->font();
    painter->drawText(-text_size.width() * 0.5, text_size.height() * 0.5, text);

    // Restore the painter's original state
    painter->restore();
}

} // namespace sketcher
} // namespace schrodinger
