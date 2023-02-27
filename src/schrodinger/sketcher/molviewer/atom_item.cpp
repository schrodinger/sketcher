#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <GraphMol/ROMol.h>

#include <Qt>
#include <QPainter>
#include <QPointF>
#include <QString>

#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a point one bond length (i.e. one unit in the RDKit coordinate
 * system) away from center, in the direction that places it furthest from any
 * of the given points
 */
QPointF best_placing_around_center(std::vector<QPointF> points, QPointF center)
{
    if (points.empty()) {
        return center + QPointF(VIEW_SCALE, 0);
    }
    std::vector<float> angles;
    // find out the angles at which each member of points is found around center
    for (auto& point : points) {
        QLineF line(center, point);
        auto angle = line.angle();
        angles.push_back(angle);
    }
    sort(angles.begin(), angles.end());
    angles.push_back(angles.front() + 360);
    // find the biggest angle interval and return a point in the middle of it
    int best_i = 0;
    auto best_angle = (angles[best_i + 1] - angles[best_i]) * 0.5;
    for (unsigned int i = 0; i < angles.size() - 1; ++i) {
        auto angle_i = (angles[i + 1] - angles[i]) * 0.5;
        if (angle_i > best_angle) {
            best_i = i;
            best_angle = angle_i;
        }
    }
    // For a single substituent, limit the angle to 120 (instead of 180)
    float max_angle = (points.size() == 1) ? 120.0 : 180.0;
    float angle_to_use = (angles[best_i + 1] - angles[best_i]) * 0.5;
    float angle_increment = std::min(max_angle, angle_to_use);
    float return_angle = angles[best_i] + angle_increment;
    QLineF line(center, center + QPointF(VIEW_SCALE, 0));
    line.setAngle(return_angle);
    return line.p2();
}

AtomItem::AtomItem(RDKit::Atom* atom, Fonts& fonts, AtomItemSettings& settings,
                   QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_atom(atom),
    m_fonts(fonts),
    m_settings(settings)
{
    setZValue(static_cast<qreal>(ZOrder::ATOM));
    m_selection_highlighting_path.addEllipse(
        QPointF(0, 0), ATOM_SELECTION_HIGHLIGHTING_RADIUS,
        ATOM_SELECTION_HIGHLIGHTING_RADIUS);
    m_predictive_highlighting_path.addEllipse(
        QPointF(0, 0), ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS,
        ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS);
    updateCachedData();
}

int AtomItem::type() const
{
    return Type;
}

void AtomItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();
    m_label_is_visible = determineLabelIsVisible();
    m_subrects.clear();
    if (m_label_is_visible) {
        m_pen.setColor(m_settings.getAtomColor(m_atom->getAtomicNum()));
        m_main_label_text = QString::fromStdString(m_atom->getSymbol());
        qreal label_width =
            m_fonts.m_main_label_fm.boundingRect(m_main_label_text).width();
        qreal label_height = m_fonts.m_main_label_fm.height() * 0.8;
        // center a label_width by label_height rectangle about (0, 0)
        m_main_label_rect = QRectF(-label_width * 0.5, -label_height * 0.5,
                                   label_width, label_height);
        m_subrects.push_back(m_main_label_rect);
    } else {
        m_main_label_text = QString();
        m_main_label_rect = QRectF();
    }

    // merge all of the subrects with the predictive highlighting path to create
    // the shape and bounding rect
    m_shape = QPainterPath(m_predictive_highlighting_path);
    for (QRectF rect : m_subrects) {
        QPainterPath rect_path;
        rect_path.addRect(rect);
        m_shape |= rect_path;
    }
    m_bounding_rect = m_shape.boundingRect();
}

bool AtomItem::determineLabelIsVisible() const
{
    if (m_settings.m_carbon_labels == CarbonLabels::ALL) {
        return true;
    }
    if (m_atom->getAtomicNum() != static_cast<unsigned int>(Element::C)) {
        return true;
    }
    // TODO: visible if isotope
    // TODO: visible if valence error
    // if (m_settings.m_show_valence_errors && ...)
    if (m_atom->getFormalCharge() != 0) {
        return true;
    }
    int num_bonds = m_atom->getDegree();
    if (num_bonds == 0) {
        return true;
    } else if (num_bonds == 1 &&
               m_settings.m_carbon_labels == CarbonLabels::TERMINAL) {
        return true;
    }
    return false;
}

const std::vector<QRectF>& AtomItem::getSubrects() const
{
    return m_subrects;
}

bool AtomItem::labelIsVisible() const
{
    return m_label_is_visible;
}

void AtomItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget)
{
    if (m_label_is_visible) {
        painter->save();
        painter->setFont(m_fonts.m_main_label_font);
        painter->setPen(m_pen);
        painter->drawText(m_main_label_rect, Qt::AlignCenter,
                          m_main_label_text);
        painter->restore();
    }
}

QPointF AtomItem::findPositionInEmptySpace(bool avoid_subrects) const
{
    QPointF lastp = scenePos();
    std::vector<QPointF> positions;
    auto& mol = m_atom->getOwningMol();
    auto& conf = mol.getConformer();
    for (const auto& neighbor : mol.atomNeighbors(m_atom)) {
        auto& pos = conf.getAtomPos(neighbor->getIdx());
        positions.push_back({pos.x, pos.y});
    }

    if (avoid_subrects) {
        for (auto rect : getSubrects()) {
            positions.push_back(rect.center() + lastp);
        }
    }
    return best_placing_around_center(positions, lastp);
}

QPointF to_scene_xy(const RDGeom::Point3D& xyz)
{
    return QPointF(xyz.x * VIEW_SCALE, -xyz.y * VIEW_SCALE);
}

} // namespace sketcher
} // namespace schrodinger
