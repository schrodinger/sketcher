#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <GraphMol/Chirality.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ROMol.h>
#include <QPainter>
#include <QPointF>
#include <QString>
#include <Qt>

#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/sketcher/molviewer/stereochemistry.h"
#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a point one bond length (i.e. one unit in the RDKit coordinate
 * system) away from the origin, in the direction that places it furthest from
 * any of the given points
 */
QPointF best_placing_around_origin(std::vector<QPointF> points)
{
    if (points.empty()) {
        return QPointF(VIEW_SCALE, 0);
    }
    QPointF origin(0.f, 0.f);
    std::vector<float> angles;
    // find out the angles at which each member of points is found around the
    // origin
    for (auto& point : points) {
        QLineF line(origin, point);
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
    QLineF line(origin, QPointF(VIEW_SCALE, 0));
    line.setAngle(return_angle);
    return line.p2();
}

/**
 * @return a bounding rect of the given label for the given font metrics, with
 * the height adjusted to better fit the actual size of the text
 */
QRect make_text_rect(QFontMetrics fm, QString label)
{
    auto rect = fm.boundingRect(label);
    return rect.adjusted(0, rect.height() * LABEL_RECT_HEIGHT_ADJUSTMENT_FACTOR,
                         0, 0);
};

AtomItem::AtomItem(const RDKit::Atom* atom, const Fonts& fonts,
                   AtomItemSettings& settings, QGraphicsItem* parent) :
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
    m_valence_error_pen = QPen(VALENCE_ERROR_BORDER_COLOR);
    m_valence_error_pen.setWidthF(VALENCE_ERROR_BORDER_WIDTH);
    m_valence_error_pen.setStyle(Qt::DotLine);
    m_valence_error_pen.setCapStyle(Qt::RoundCap);
    m_valence_error_brush = QBrush(VALENCE_ERROR_AREA_COLOR);

    // setup chirality pen
    m_chirality_pen = QPen(CHIRALITY_LABEL_COLOR);

    updateCachedData();
}

int AtomItem::type() const
{
    return Type;
}

HsDirection AtomItem::findHsDirection() const
{
    auto& mol = m_atom->getOwningMol();
    auto& conf = mol.getConformer();
    switch (m_atom->getDegree()) {
        case 0: {
            // if there's no bonds, draw Hs to the left for chalcogens and
            // halogens and to the right for any other atom
            const auto* table = RDKit::PeriodicTable::getTable();
            auto n_electrons = table->getNouterElecs(m_atom->getAtomicNum());
            return ((n_electrons == 6 || n_electrons == 7)
                        ? HsDirection::LEFT
                        : HsDirection::RIGHT);
        }
        case 1: {
            // with 1 bond, draw the Hs horizontally, on the opposite side of
            // the bound atom
            auto neighbor = mol.atomNeighbors(m_atom).begin();
            auto pos = conf.getAtomPos((*neighbor)->getIdx()) -
                       conf.getAtomPos(m_atom->getIdx());
            return (pos.x > 0 ? HsDirection::LEFT : HsDirection::RIGHT);
        }
        default: {
            // with two or more bonds, draw the Hs vertically or horizontally,
            // as far away as possible from bonded atoms
            auto position = findPositionInEmptySpace(false);

            if (std::abs(position.x()) > std::abs(position.y())) {
                return (position.x() > 0) ? HsDirection::RIGHT
                                          : HsDirection::LEFT;
            } else {
                return (position.y() > 0) ? HsDirection::DOWN : HsDirection::UP;
            }
        }
    }
}

void AtomItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();
    clearLabels();

    m_label_is_visible = determineLabelIsVisible();
    if (m_label_is_visible) {
        m_pen.setColor(m_settings.getAtomColor(m_atom->getAtomicNum()));
        m_main_label_text = QString::fromStdString(m_atom->getSymbol());

        // if there's a user-set label, override m_main_label_text
        if (!m_user_label.isEmpty()) {
            m_main_label_text = m_user_label;
        }

        m_main_label_rect =
            make_text_rect(m_fonts.m_main_label_fm, m_main_label_text);
        m_main_label_rect.moveCenter(QPointF(0, 0));
        if (m_user_label.isEmpty()) {
            updateIsotopeLabel();
            updateChargeAndRadicalLabel();
            updateHsLabel();
            positionLabels();
        }

        for (auto rect : getLabelRects()) {
            if (rect.isValid()) {
                m_subrects.push_back(rect);
            }
        }
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
    updateChiralityLabel();
}

void AtomItem::updateChiralityLabel()
{
    auto chirality = get_atom_chirality_label(*m_atom);
    if (!chirality.isEmpty()) {
        // add parentheses around chiral label
        chirality = "(" + chirality + ")";
        // TODO: enhanced stereo labels - SKETCH-1960
    }
    m_chirality_label_text = chirality;
    auto chirality_label_position =
        findPositionInEmptySpace(true) * CHIRALITY_LABEL_DISTANCE_RATIO;
    m_chirality_label_rect =
        make_text_rect(m_fonts.m_chirality_fm, m_chirality_label_text);
    m_chirality_label_rect.moveCenter(chirality_label_position);
}

void AtomItem::clearLabels()
{
    for (auto string : {m_main_label_text, m_charge_and_radical_label_text,
                        m_H_count_label_text, m_chirality_label_text}) {
        string = QString();
    }
    for (auto rect : getLabelRects()) {
        rect = QRectF();
    }
    m_subrects.clear();
}

void AtomItem::updateIsotopeLabel()
{
    auto isotope = m_atom->getIsotope();
    // when no isotope is set RDKit returns 0
    if (isotope != 0) {
        m_isotope_label_text = QString::number(isotope);
        m_isotope_rect =
            make_text_rect(m_fonts.m_subscript_fm, m_isotope_label_text);
    }
}

void AtomItem::updateChargeAndRadicalLabel()
{
    QString radical = "•"; // U+2022
    auto radical_electrons_count = m_atom->getNumRadicalElectrons();
    auto charge = m_atom->getFormalCharge();
    QString radical_label("");
    if (radical_electrons_count > 0) {
        auto parentheses = radical_electrons_count > 1 && charge != 0;
        radical_label = (parentheses ? "(" : "") +
                        (radical_electrons_count > 1
                             ? QString::number(radical_electrons_count)
                             : "") +
                        radical + (parentheses ? ")" : "");
    }
    QString sign("+");
    QString charge_label("");

    if (charge < 0) {
        charge = -charge;
        // en dash
        sign = QString("–");
    }
    if (charge == 1) {
        charge_label = sign;
    } else if (charge > 1) {
        charge_label = QString::number(charge) + sign;
    }
    auto text_label =
        radical_label +
        (radical_label.isEmpty() || charge_label.isEmpty() ? "" : " ") +
        charge_label;
    if (!text_label.isEmpty()) {
        m_charge_and_radical_label_text = text_label;
        m_charge_and_radical_rect = make_text_rect(
            m_fonts.m_subscript_fm, m_charge_and_radical_label_text);
    }
}

void AtomItem::updateHsLabel()
{
    bool includeNeighbors = false;
    auto H_count = m_atom->getTotalNumHs(includeNeighbors);

    if (H_count > 0) {
        m_H_label_rect = make_text_rect(m_fonts.m_main_label_fm, "H");
    }
    if (H_count > 1) {
        m_H_count_label_text = QString::number(H_count);
        m_H_count_label_rect =
            make_text_rect(m_fonts.m_subscript_fm, m_H_count_label_text);
    }
}

void AtomItem::positionLabels()
{
    auto label_rect = m_main_label_rect;
    m_isotope_rect.translate(
        QPointF(m_main_label_rect.left(), m_main_label_rect.center().y()) -
        m_isotope_rect.bottomRight());
    m_charge_and_radical_rect.translate(
        QPointF(m_main_label_rect.right(), m_main_label_rect.center().y()) -
        m_charge_and_radical_rect.bottomLeft());
    label_rect = label_rect.united(m_isotope_rect);
    auto H_direction = findHsDirection();
    // charge and radical text is next to the main label, on the right, unless
    // hydrogens are drawn to the right, in which case they come before charge
    // and radicals. (e.g. H3O+ but NH4+)
    if (H_direction != HsDirection::RIGHT) {
        label_rect = label_rect.united(m_charge_and_radical_rect);
    }
    auto SPACER = m_main_label_rect.width() * LABEL_SPACER_RATIO;
    if (!m_H_count_label_rect.isValid() && H_direction == HsDirection::LEFT) {
        label_rect.adjust(-SPACER, 0, 0, 0);
    }
    if (m_H_label_rect.isValid() && H_direction == HsDirection::RIGHT) {
        label_rect.adjust(0, 0, SPACER, 0);
    }
    switch (H_direction) {
        case HsDirection::LEFT:
            m_H_count_label_rect.translate(
                label_rect.bottomLeft() -
                QPointF(m_H_count_label_rect.right(),
                        m_H_count_label_rect.center().y()));
            label_rect.setLeft(m_H_count_label_rect.left());
            m_H_label_rect.translate(label_rect.bottomLeft() -
                                     m_H_label_rect.bottomRight());
            break;

        case HsDirection::RIGHT:
            m_H_label_rect.translate(label_rect.bottomRight() -
                                     m_H_label_rect.bottomLeft());
            label_rect.setRight(m_H_label_rect.right());
            m_H_count_label_rect.translate(
                label_rect.bottomRight() -
                QPointF(m_H_count_label_rect.left(),
                        m_H_count_label_rect.center().y()));
            label_rect.setRight(m_H_count_label_rect.right());
            m_charge_and_radical_rect.translate(
                m_H_count_label_rect.right() - m_charge_and_radical_rect.left(),
                0);

            break;
        case HsDirection::UP: {
            auto label_and_H_height =
                m_main_label_rect.height() + m_H_label_rect.height();
            m_H_label_rect.moveCenter(m_main_label_rect.center() +
                                      QPointF(0, label_and_H_height * 0.5));
            m_H_count_label_rect.translate(
                m_H_label_rect.bottomRight() -
                QPointF(m_H_count_label_rect.left(),
                        m_H_count_label_rect.center().y()));
            QPointF offset(0, label_rect.top() - m_H_count_label_rect.bottom());
            if (m_H_count_label_rect.isValid() &&
                m_charge_and_radical_rect.isValid()) {
                offset -= QPoint(0, SPACER);
            }
            m_H_label_rect.translate(offset);
            m_H_count_label_rect.translate(offset);
            break;
        }
        case HsDirection::DOWN: {
            auto label_and_H_height =
                m_main_label_rect.height() + m_H_label_rect.height();
            m_H_label_rect.moveCenter(m_main_label_rect.center() +
                                      QPointF(0, label_and_H_height * 0.5));
            m_H_count_label_rect.translate(
                m_H_label_rect.bottomRight() -
                QPointF(m_H_count_label_rect.left(),
                        m_H_count_label_rect.center().y()));
            break;
        }
    }
}

bool AtomItem::determineLabelIsVisible() const
{
    // note that this method does not take m_valence_error_is_visible into
    // account.  If that value is true (and up-to-date), then the label should
    // be visible.
    if (!m_user_label.isEmpty()) {
        return true;
    }
    if (m_settings.m_carbon_labels == CarbonLabels::ALL) {
        return true;
    }
    if (m_atom->getAtomicNum() != static_cast<unsigned int>(Element::C)) {
        return true;
    }
    if (m_atom->getIsotope() != 0) {
        return true;
    }
    if (m_atom->getNumRadicalElectrons()) {
        return true;
    }
    if (m_atom->getFormalCharge() != 0) {
        return true;
    }
    int num_bonds = m_atom->getDegree();
    if (num_bonds == 0) {
        return true;
    } else if (num_bonds == 1 &&
               m_settings.m_carbon_labels == CarbonLabels::TERMINAL) {
        return true;
    } else if (num_bonds == 2) {
        // show label of carbons with two double bonds
        auto bonds = m_atom->getOwningMol().atomBonds(m_atom);
        auto is_double = [](auto bond) {
            return bond->getBondType() == RDKit::Bond::DOUBLE;
        };
        if (std::all_of(bonds.begin(), bonds.end(), is_double)) {
            return true;
        }
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
    if (m_valence_error_is_visible) {
        painter->save();
        painter->setPen(m_valence_error_pen);
        painter->setBrush(m_valence_error_brush);
        painter->drawEllipse(m_main_label_rect.adjusted(
            -VALENCE_ERROR_AREA_BORDER, -VALENCE_ERROR_AREA_BORDER,
            VALENCE_ERROR_AREA_BORDER, VALENCE_ERROR_AREA_BORDER));
        painter->restore();
    }
    if (m_label_is_visible) {
        painter->save();
        painter->setFont(m_fonts.m_main_label_font);
        painter->setPen(m_pen);
        painter->drawText(m_main_label_rect, Qt::AlignCenter,
                          m_main_label_text);
        // drawText only paints when the rect is valid
        painter->drawText(m_H_label_rect, Qt::AlignCenter, "H");
        painter->setFont(m_fonts.m_subscript_font);
        painter->drawText(m_H_count_label_rect, Qt::AlignCenter,
                          m_H_count_label_text);
        painter->drawText(m_charge_and_radical_rect, Qt::AlignCenter,
                          m_charge_and_radical_label_text);
        painter->drawText(m_isotope_rect, Qt::AlignCenter,
                          m_isotope_label_text);
        painter->restore();
    }
    painter->save();
    painter->setPen(m_chirality_pen);
    painter->setFont(m_fonts.m_chirality_font);
    painter->drawText(m_chirality_label_rect, m_chirality_label_text);
    painter->restore();
}

const RDKit::Atom* AtomItem::getAtom() const
{
    return m_atom;
}

QPointF AtomItem::findPositionInEmptySpace(bool avoid_subrects) const
{

    auto& mol = m_atom->getOwningMol();
    auto& conf = mol.getConformer();
    std::vector<QPointF> positions;
    auto& this_pos = conf.getAtomPos(m_atom->getIdx());
    QPointF lastp = to_scene_xy(this_pos);

    for (const auto& neighbor : mol.atomNeighbors(m_atom)) {
        auto& pos = conf.getAtomPos(neighbor->getIdx());
        positions.push_back(to_scene_xy(pos) - lastp);
    }

    if (avoid_subrects) {
        for (auto rect : getSubrects()) {
            if (rect == m_main_label_rect)
                continue;
            positions.push_back(rect.topLeft());
            positions.push_back(rect.topRight());
            positions.push_back(rect.bottomLeft());
            positions.push_back(rect.bottomRight());
            positions.push_back(rect.center());
        }
    }
    return best_placing_around_origin(positions);
}

std::vector<QRectF> AtomItem::getLabelRects() const
{
    return {m_main_label_rect,    m_isotope_rect, m_charge_and_radical_rect,
            m_H_count_label_rect, m_H_label_rect, m_chirality_label_rect};
}

QPointF to_scene_xy(const RDGeom::Point3D& xyz)
{
    return QPointF(xyz.x * VIEW_SCALE, -xyz.y * VIEW_SCALE);
}

bool AtomItem::determineValenceErrorIsVisible() const
{
    if (!m_settings.m_valence_errors_shown) {
        return false;
    }
    auto atomic_number = m_atom->getAtomicNum();
    const auto* table = RDKit::PeriodicTable::getTable();

    /*  TODO: SKETCH-1917 non-elements, or elements signifying any valence is
       permitted, or attached to any query bond if (isQuery() ||
       table->getDefaultValence(atomic_number) == -1 ||
           std::any_of(m_bonds.begin(), m_bonds.end(),
                       [](const auto& b) { return b->isQuery(); })) {
           return false;
       }
     */

    // There is a valence error if the current value is not permitted
    auto allowed_valence = table->getValenceList(atomic_number);

    // RDKit's getTotalValence() doesn't consider charges or unpaired electrons
    // -- it returns the total bond order instead of the total valence.
    auto current_valence = m_atom->getTotalValence() -
                           m_atom->getFormalCharge() +
                           m_atom->getNumRadicalElectrons();
    return std::find(allowed_valence.begin(), allowed_valence.end(),
                     current_valence) == allowed_valence.end();
}

} // namespace sketcher
} // namespace schrodinger
