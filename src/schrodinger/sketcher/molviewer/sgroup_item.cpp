#include "schrodinger/sketcher/rdkit/sgroup.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/rdkit/s_group_constants.h"
#include <GraphMol/ROMol.h>
#include <QMarginsF>
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
    show();
}

std::pair<QPointF, bool> SGroupItem::getFieldDataDisplayInfo() const
{
    std::string fieldDisp;
    if (m_sgroup.getPropIfPresent("FIELDDISP", fieldDisp)) {
        QString displacement = QString::fromStdString(fieldDisp);
        auto disp_x = displacement.sliced(0, 10).toFloat();
        auto disp_y = displacement.sliced(10, 10).toFloat();
        // if the string was not a valid float, the toFloat() function will
        // return 0, which is what we want, so no need to check for errors
        return {QPointF(disp_x, disp_y), displacement[26] == 'R'};
    }
    return {QPointF(), false};
}

/**
 * @return The repeat pattern label text for the SGroup, if any.
 * NOTE: Users expect head-to-tail (ht) to not be rendered.
 */
static QString get_repeat_connect_text(const RDKit::SubstanceGroup& sgroup)
{
    QString text;
    auto repeat_str = get_repeat_pattern_label(sgroup);
    if (REPEATPATTERN_TO_RDKITSTRING_BIMAP.right.count(repeat_str)) {
        auto repeat = REPEATPATTERN_TO_RDKITSTRING_BIMAP.right.at(repeat_str);
        if (repeat != RepeatPattern::HEAD_TO_TAIL) {
            text = QString::fromStdString(repeat_str).toLower();
        }
    }
    return text;
}

/**
 * @return The sgroup type label text for the SGroup, if any.
 * NOTE: Accounts for defaults for both SRU and COP sgroup types.
 */
static QString get_repeat_label_text(const RDKit::SubstanceGroup& sgroup)
{
    QString text;
    auto type_str = get_sgroup_type(sgroup);
    if (SUBGROUPTYPE_TO_RDKITSTRING_BIMAP.right.count(type_str)) {
        auto type = SUBGROUPTYPE_TO_RDKITSTRING_BIMAP.right.at(type_str);
        auto label = get_polymer_label(sgroup);
        if (type == SubgroupType::SRU_POLYMER && label.empty()) {
            text = "n";
        } else if (type == SubgroupType::COPOLYMER && label.empty()) {
            text = "co";
        } else {
            text = QString::fromStdString(label);
        }
    }
    return text;
}

/**
 * @return The FIELDDATA text for the SGroup, if any. This text is displayed
 * next to the SGroup atoms.
 */
static QString get_field_data_text(const RDKit::SubstanceGroup& sgroup)
{
    QString text;
    std::string typ;
    if (sgroup.getPropIfPresent("TYPE", typ) && typ == "DAT") {
        if (sgroup.hasProp("DATAFIELDS")) {
            auto dfs = sgroup.getProp<std::vector<std::string>>("DATAFIELDS");
            for (const auto& df : dfs) {
                text += df + "|";
            }
            text.chop(1);
        }
    }
    return text;
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
    if (!m_field_data_text.isEmpty()) {
        painter->drawText(m_field_data_text_rect.bottomLeft(),
                          m_field_data_text);
    }
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

void SGroupItem::updateCachedData()
{
    prepareGeometryChange();
    m_brackets_path = getBracketPath();

    m_repeat = get_repeat_connect_text(m_sgroup);
    m_label = get_repeat_label_text(m_sgroup);
    m_field_data_text = get_field_data_text(m_sgroup);

    // until https://github.com/rdkit/rdkit/issues/7829 is resolved, we need to
    // place the field data text manually
    auto [field_data_coords, field_data_coords_are_relative] =
        getFieldDataDisplayInfo();
    m_field_data_text_rect =
        m_fonts.m_sgroup_fm.tightBoundingRect(m_field_data_text);
    // for now we ignore the relative flag and always place the field data text
    // relative to the first atom in the sgroup
    Q_UNUSED(field_data_coords_are_relative);
    auto center = m_sgroup.getOwningMol().getConformer().getAtomPos(
        m_sgroup.getAtoms()[0]);
    m_field_data_text_rect.moveCenter(to_scene_xy(center) + field_data_coords);
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

    // generate the shape and highlights after everything else, since they
    // incorporate the above values
    m_shape = getShapeWithWidth(S_GROUP_HIGHLIGHT_PADDING);
    m_bounding_rect = m_shape.boundingRect();
    m_selection_highlighting_path = m_shape;

    auto atom_and_bond_pred_path =
        get_predictive_highlighting_path_for_s_group_atoms_and_bonds(m_sgroup);
    m_predictive_highlighting_path =
        m_shape + mapFromScene(atom_and_bond_pred_path);
}

QPainterPath SGroupItem::getShapeWithWidth(qreal width) const
{
    QPainterPathStroker stroker;
    stroker.setWidth(width);
    auto path = stroker.createStroke(m_brackets_path);
    // we intentionally add a little extra highlighting to the right of the
    // labels.  Otherwise, it looks unbalanced because the brackets are to the
    // left of the labels.
    const QMarginsF expand_labels_by(width / 2.0, width / 2.0, width,
                                     width / 2.0);
    if (!m_label.isEmpty()) {
        auto expanded_label =
            m_untransformed_label_rect.marginsAdded(expand_labels_by);
        auto transformed_label = m_labels_transform.map(expanded_label);
        path.addPolygon(transformed_label);
    }
    if (!m_repeat.isEmpty()) {
        auto expanded_repeat =
            m_untransformed_repeat_rect.marginsAdded(expand_labels_by);
        auto transformed_repeat = m_labels_transform.map(expanded_repeat);
        path.addPolygon(transformed_repeat);
    }
    if (!m_field_data_text.isEmpty()) {
        path.addRect(m_field_data_text_rect);
    }
    return path;
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
    // DAT substance groups on CG mols have no brackets since they apply to the
    // entire structure
    if (brackets.size() == 0) {
        return {QLineF(), QPointF()};
    }

    typedef std::array<RDGeom::Point3D, 3> Bracket;
    auto rightmost_x = [](const Bracket& bracket) {
        return std::max(bracket[0].x, bracket[1].x);
    };
    auto rightmost_bracket = *(std::max_element(
        brackets.begin(), brackets.end(),
        [rightmost_x](const Bracket& bracket_1, const Bracket& bracket_2) {
            return rightmost_x(bracket_1) < rightmost_x(bracket_2);
        }));
    // return a QLineF with both points of the brackets in scene coordinates.
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