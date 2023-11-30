#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/rdkit_extensions/sgroup.h"
#include <GraphMol/ROMol.h>
#include <QPainter>

namespace schrodinger
{
namespace sketcher
{

SGroupItem::SGroupItem(const RDKit::SubstanceGroup& sgroup, const Fonts& fonts,
                       QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_sgroup(sgroup),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::SGROUP));
    updateCachedData();
}

const RDKit::SubstanceGroup* SGroupItem::getSubstanceGroup() const
{
    return &m_sgroup;
}

void SGroupItem::paint(QPainter* painter,
                       const QStyleOptionGraphicsItem* option, QWidget* widget)
{
    Q_UNUSED(option);
    Q_UNUSED(widget);
    painter->save();
    painter->setPen(Qt::black);
    painter->drawPath(m_brackets_path);
    painter->setFont(m_fonts.m_sgroup_font);
    painter->setTransform(m_labels_transform * painter->transform());
    painter->drawText(m_untransformed_label_rect.bottomLeft(), m_label);
    painter->drawText(m_untransformed_repeat_rect.bottomLeft(), m_repeat);
    painter->restore();
}

int SGroupItem::type() const
{
    return Type;
}

QPainterPath SGroupItem::computePredictiveHighlightingPath() const
{
    QPainterPath path = shape();
    auto& molecule = m_sgroup.getOwningMol();
    auto scene_ptr = dynamic_cast<Scene*>(scene());
    if (scene_ptr) {
        for (auto atom_idx : m_sgroup.getAtoms()) {
            auto atom_item = scene_ptr->getAtomItemForAtom(
                molecule.getAtomWithIdx(atom_idx));
            if (atom_item) {
                path.addPath(mapFromItem(
                    atom_item, atom_item->predictiveHighlightingPath()));
            }
        }
    }
    return path;
}

void SGroupItem::updateCachedData()
{
    prepareGeometryChange();
    m_repeat = QString::fromStdString(
        rdkit_extensions::get_repeat_pattern_label(m_sgroup));
    m_label =
        QString::fromStdString(rdkit_extensions::get_polymer_label(m_sgroup));
    m_brackets_path = getBracketPath();
    m_predictive_highlighting_path = computePredictiveHighlightingPath();
    auto [positions, displacement] = getPositionsForLabels();
    auto brackets_half_height = positions.length() * 0.5;
    auto translation_offset = (positions.p1() + positions.p2()) * 0.5;
    m_untransformed_label_rect = m_fonts.m_sgroup_fm.tightBoundingRect(m_label);
    m_untransformed_repeat_rect =
        m_fonts.m_sgroup_fm.tightBoundingRect(m_repeat);
    qreal label_x_offset = LABEL_DISTANCE_FROM_BRACKETS * VIEW_SCALE;
    qreal label_y_offset = brackets_half_height;
    QPointF label_rect_untransformed_offset(label_x_offset, label_y_offset);
    QPointF repeat_rect_untransformed_offset(label_x_offset, -label_y_offset);
    // move the label rects so the middle of their left side is label_x_offset
    // away from the brackets ends
    m_untransformed_label_rect.translate(
        label_rect_untransformed_offset -
        (m_untransformed_label_rect.bottomLeft() +
         m_untransformed_label_rect.topLeft()) *
            0.5);
    m_untransformed_repeat_rect.translate(
        repeat_rect_untransformed_offset -
        (m_untransformed_repeat_rect.bottomLeft() +
         m_untransformed_repeat_rect.topLeft()) *
            0.5);
    auto rotation_angle =
        QLineF(QPointF(0, 0), displacement).angleTo(QLineF(0, 0, 1, 0));
    QTransform transform;
    transform.translate(translation_offset.x(), translation_offset.y());
    transform.rotate(rotation_angle);
    m_labels_transform = transform;
}

QPainterPath SGroupItem::shape() const
{
    QPainterPathStroker stroker;
    stroker.setWidth(4);
    auto path = stroker.createStroke(m_brackets_path);
    path.addPolygon(m_labels_transform.map(m_untransformed_label_rect));
    path.addPolygon(m_labels_transform.map(m_untransformed_repeat_rect));
    return path;
}

QRectF SGroupItem::boundingRect() const
{
    return shape().boundingRect();
}

QPainterPath SGroupItem::getBracketPath() const
{
    QPainterPath path;
    for (auto bracket : m_sgroup.getBrackets()) {
        auto short_side = computeShortSideForBracket(bracket[0], bracket[1]);
        path.moveTo(to_scene_xy(bracket[0] + short_side));
        path.lineTo(to_scene_xy(bracket[0]));
        path.lineTo(to_scene_xy(bracket[1]));
        path.lineTo(to_scene_xy(bracket[1] + short_side));
    }
    return path;
}

std::pair<QLineF, QPointF> SGroupItem::getPositionsForLabels() const
{
    // find the bracket that has the point with the highest x
    auto brackets = m_sgroup.getBrackets();
    typedef std::array<RDGeom::Point3D, 3> Bracket;
    auto rightmost_x = [](const Bracket& bracket) {
        return std::max(bracket[0].x, bracket[1].x);
    };
    auto rightmost_bracket = *(std::max_element(
        brackets.begin(), brackets.end(),
        [rightmost_x](const Bracket& bracket_1, const Bracket& bracket_2) {
            return rightmost_x(bracket_1) < rightmost_x(bracket_2);
        }));
    // return a QLineF with both points of the brackets in scene
    // coordinates.
    QLineF out(to_scene_xy(rightmost_bracket[0]),
               to_scene_xy(rightmost_bracket[1]));

    // Make sure that the highest y point is the first one
    if (out.dy() > 0) {
        out = QLineF(out.p2(), out.p1());
    }
    QPointF direction = to_scene_xy(-computeShortSideForBracket(
        rightmost_bracket[0], rightmost_bracket[1]));
    return {out, direction};
}

RDGeom::Point3D
SGroupItem::computeShortSideForBracket(const RDGeom::Point3D& point1,
                                       const RDGeom::Point3D& point2) const
{
    QLineF bracket_line(point1.x, point1.y, point2.x, point2.y);
    auto molecule = m_sgroup.getOwningMol();
    auto connection_bonds = m_sgroup.getBonds();
    auto atoms = m_sgroup.getAtoms();
    for (auto idx : connection_bonds) {
        auto connection_bond = molecule.getBondWithIdx(idx);
        auto atom1 = connection_bond->getBeginAtomIdx();
        auto atom2 = connection_bond->getEndAtomIdx();
        auto pos1 = molecule.getConformer().getAtomPos(atom1);
        auto pos2 = molecule.getConformer().getAtomPos(atom2);
        QLineF bond_line(pos1.x, pos1.y, pos2.x, pos2.y);
        QPointF intersection_point;
        if (bracket_line.intersects(bond_line, &intersection_point)) {
            auto atom_pos =
                (find(atoms.begin(), atoms.end(), atom1) != atoms.end() ? pos1
                                                                        : pos2);
            auto bracket_dir = point1 - point2;
            auto normal = RDGeom::Point3D(-bracket_dir.y, bracket_dir.x, 0.0);
            normal.normalize();
            if (normal.dotProduct(atom_pos -
                                  RDGeom::Point3D(intersection_point.x(),
                                                  intersection_point.y(), 0)) <
                0) {
                normal *= -1;
            }
            return normal * BRACKETS_SHORT_SIDE;
        }
    }
    return RDGeom::Point3D();
}

} // namespace sketcher
} // namespace schrodinger