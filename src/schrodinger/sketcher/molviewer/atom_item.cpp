#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/RDGeneral/types.h>

#include <QPainter>
#include <QPointF>
#include <QString>
#include <Qt>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"
#include "schrodinger/rdkit_extensions/variable_attachment_bond.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/queries.h"

namespace schrodinger
{
namespace sketcher
{

bool has_enhanced_query_properties(const RDKit::Atom* atom)
{
    if (!atom->hasQuery()) {
        return false;
    }
    auto no_advanced_props = AtomQueryProperties();
    auto props = read_properties_from_atom(atom);
    auto query_props = *static_cast<const AtomQueryProperties*>(props.get());
    no_advanced_props.query_type = query_props.query_type;
    no_advanced_props.allowed_list = query_props.allowed_list;
    no_advanced_props.charge = query_props.charge;
    no_advanced_props.unpaired_electrons = query_props.unpaired_electrons;
    no_advanced_props.isotope = query_props.isotope;
    no_advanced_props.enhanced_stereo = query_props.enhanced_stereo;
    no_advanced_props.element = query_props.element;

    return query_props != no_advanced_props;
}

/**
 * @return a bounding rect of the given label for the given font metrics
 */
QRectF make_text_rect(QFontMetricsF fm, QString label)

{
    // if the label is empty return an invalid QRectF
    if (label.isEmpty()) {
        return QRectF();
    }

    auto rect = fm.tightBoundingRect(label);
    // add some margin around the text
    return rect.adjusted(-TEXT_RECT_MARGIN, -TEXT_RECT_MARGIN, TEXT_RECT_MARGIN,
                         TEXT_RECT_MARGIN);
};

AtomItem::AtomItem(const RDKit::Atom* atom, const Fonts& fonts,
                   const AtomDisplaySettings& settings, QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_atom(atom),
    m_fonts(fonts),
    m_settings(settings)
{
    setZValue(static_cast<qreal>(ZOrder::ATOM));
    m_valence_error_pen = QPen(VALENCE_ERROR_BORDER_COLOR);
    m_valence_error_pen.setWidthF(VALENCE_ERROR_BORDER_WIDTH);
    m_valence_error_pen.setStyle(Qt::DotLine);
    m_valence_error_pen.setCapStyle(Qt::RoundCap);
    m_valence_error_brush = QBrush(VALENCE_ERROR_AREA_COLOR);

    // setup chirality pen
    m_chirality_pen = QPen(CHIRALITY_LABEL_COLOR);

    // use the default atom color for the squiggle pen (used for attachment
    // points)
    m_squiggle_pen = QPen(settings.getAtomColor(-1));

    m_element_list_line_pen = QPen(ELEMENT_LIST_LINE_COLOR);
    m_squiggle_pen.setWidthF(BOND_DEFAULT_PEN_WIDTH);
    m_squiggle_pen.setCapStyle(Qt::RoundCap);

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

void AtomItem::setHideStereoLabels(bool hide)
{
    m_hide_stereo_labels = hide;
}

void AtomItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();
    clearLabels();
    m_selection_highlighting_path =
        get_selection_highlighting_path_for_atom(m_atom);
    m_predictive_highlighting_path =
        get_predictive_highlighting_path_for_atom(m_atom);
    m_squiggle_pen.setWidthF(m_settings.m_squiggle_pen_width);

    bool needs_additional_labels, show_Hs;
    std::tie(m_main_label_text, m_squiggle_path, m_label_is_visible,
             m_valence_error_is_visible, needs_additional_labels, show_Hs,
             m_query_label_text) = determineLabelType();
    if (m_label_is_visible) {
        m_pen.setColor(m_settings.getAtomColor(m_atom->getAtomicNum()));
        m_main_label_rect =
            make_text_rect(m_fonts.m_main_label_fm, m_main_label_text);
        m_main_label_rect.moveCenter(QPointF(0, 0));
        if (needs_additional_labels) {
            updateIsotopeLabel();
            updateChargeAndRadicalLabel();
            if (show_Hs) {
                updateHsLabel();
            }
            updateMappingLabel();
            positionLabels();
        }
    }
    if (!m_query_label_text.isEmpty()) {
        m_query_label_rect =
            make_text_rect(m_fonts.m_element_list_fm, m_query_label_text);
        // find a position for the element list label
        auto position =
            findPositionInEmptySpace(false) * ELEMENT_LIST_LABEL_DISTANCE_RATIO;
        m_query_label_rect.moveCenter(position);
    }
    if (m_settings.m_stereo_labels_shown && !m_hide_stereo_labels) {
        updateChiralityLabel();
    }
    for (auto rect : getLabelRects()) {
        if (rect.isValid()) {
            m_subrects.push_back(rect);
        }
    }
    // merge all of the subrects with the predictive highlighting path to create
    // the shape and bounding rect
    m_shape = QPainterPath(m_predictive_highlighting_path);
    m_shape.setFillRule(Qt::WindingFill);
    for (QRectF rect : m_subrects) {
        m_shape.addRect(rect);
    }
    m_bounding_rect = m_shape.boundingRect();
}

QString AtomItem::getQueryLabel() const
{
    QString query_text;
    auto props = read_properties_from_atom(m_atom);
    auto query_props = *static_cast<const AtomQueryProperties*>(props.get());

    // if this is an allowed or disallowed list with no added features, we can
    // use a simplified label, otherwise we need to show the SMARTS
    auto no_advanced_props = AtomQueryProperties();
    no_advanced_props.query_type = query_props.query_type;
    no_advanced_props.allowed_list = query_props.allowed_list;
    bool has_advanced_props = query_props != no_advanced_props;

    if ((query_props.query_type == QueryType::ALLOWED_LIST ||
         query_props.query_type == QueryType::NOT_ALLOWED_LIST) &&
        !has_advanced_props) {
        auto list = join_all_atomic_symbols(query_props.allowed_list);
        if (!list.isEmpty()) {
            if (query_props.query_type == QueryType::NOT_ALLOWED_LIST) {
                list = "!" + list;
            }
            query_text = QString("[%1]").arg(list);
        }
    } else {
        auto smarts = (query_props.query_type == QueryType::SMARTS
                           ? query_props.smarts_query
                           : get_atom_smarts(m_atom));
        query_text = QString("\"%1\"").arg(QString::fromStdString(smarts));

        // if the string is too long, truncate it and add an ellipsis
        query_text = m_fonts.m_element_list_fm.elidedText(
            query_text, Qt::ElideRight, MAX_SMARTS_LABEL_LENGTH);
    }
    return query_text;
}
std::tuple<QString, QPainterPath, bool, bool, bool, bool, QString>
AtomItem::determineLabelType() const
{
    std::string main_label_text = "";
    QPainterPath squiggle_path;
    bool valence_error_is_visible = false;
    bool label_is_visible = true;
    bool needs_additional_labels = false;
    bool H_labels_are_visible = true;
    QString query_label_text;

    // if there's a user-set display label, always use that.
    if (m_atom->hasProp(RDKit::common_properties::_displayLabel)) {
        main_label_text = m_atom->getProp<std::string>(
            RDKit::common_properties::_displayLabel);
    } else if (m_atom->getAtomicNum() !=
                   rdkit_extensions::DUMMY_ATOMIC_NUMBER &&
               m_atom->hasProp(RDKit::common_properties::atomLabel)) {

        // if there's no user-set, but an atomLabel is present, always
        // display that on non-dummy atoms. Query atoms are dealt with below
        main_label_text =
            m_atom->getProp<std::string>(RDKit::common_properties::atomLabel);

    } else if (m_atom->hasQuery() && !m_atom->getQueryType().empty()) {
        // the query type is set for wildcards but not most other queries
        main_label_text = m_atom->getQueryType();
    } else if (m_atom->getAtomicNum() ==
               rdkit_extensions::DUMMY_ATOMIC_NUMBER) {
        if (is_attachment_point(m_atom)) {
            label_is_visible = false;
            squiggle_path = getWavyLine();
        } else if (auto r_group_num =
                       rdkit_extensions::get_r_group_number(m_atom)) {

            main_label_text = "R" + std::to_string(r_group_num.value());
        } else if (m_atom->hasProp(RDKit::common_properties::atomLabel)) {

            // display the atomLabel if present. Note that this option needs to
            // go after R-groups and attachment points because they both use the
            // atomLabel property but we don't want to display it for them
            main_label_text = m_atom->getProp<std::string>(
                RDKit::common_properties::atomLabel);
        } else if (rdkit_extensions::is_dummy_atom_for_variable_attachment_bond(
                       m_atom)) {

            // dummy atoms for variable attachment bonds aren't shown, but we
            // still need to set label_is_visible to false here so that the
            // bound BondItem knows not to make space for the "*" label
            label_is_visible = false;
        } else if (m_atom->hasQuery()) {
            // this needs to be after R groups and attachment points, since it
            // would trigger for both but we need to special-case them
            label_is_visible = false;
            auto props = read_properties_from_atom(m_atom);
            if (props->isQuery()) {
                query_label_text = getQueryLabel();
            }
        } else {
            // Unrecognized dummy atom.  Display any user-set labels
            main_label_text = m_atom->getSymbol();
        }
    } else if (m_atom->hasQuery()) {

        main_label_text = m_atom->getSymbol();
        needs_additional_labels = true;
        // hide implicit Hs on queries
        H_labels_are_visible = false;

        // if advanced properties are present on the query, show the SMARTS
        if (has_enhanced_query_properties(m_atom)) {
            label_is_visible = false;
            needs_additional_labels = false;
            auto props = read_properties_from_atom(m_atom);
            if (props->isQuery()) {
                query_label_text = getQueryLabel();
            }
        }
    } else {
        // a "normal" atom, i.e. an atom that represents an element
        valence_error_is_visible = determineValenceErrorIsVisible();
        label_is_visible =
            valence_error_is_visible || determineLabelIsVisible();
        if (label_is_visible) {
            main_label_text = m_atom->getSymbol();
            needs_additional_labels = true;
        }
    }
    // if the atom has no bonds and we are showing the element list or smarts
    // label, show an asterisk as the main label
    if (!query_label_text.isEmpty() && m_atom->getDegree() == 0 &&
        !label_is_visible) {
        main_label_text = "*";
        label_is_visible = true;
    }

    return {QString::fromStdString(main_label_text),
            squiggle_path,
            label_is_visible,
            valence_error_is_visible,
            needs_additional_labels,
            H_labels_are_visible,
            query_label_text};
}

QPainterPath AtomItem::getWavyLine() const
{
    qreal angle = get_attachment_point_line_angle(m_atom);
    return get_wavy_line_path(ATTACHMENT_POINT_SQUIGGLE_NUMBER_OF_WAVES,
                              ATTACHMENT_POINT_SQUIGGLE_WIDTH_PER_WAVE,
                              ATTACHMENT_POINT_SQUIGGLE_HEIGHT, angle);
}

void AtomItem::updateChiralityLabel()
{
    auto strip_abs = !m_settings.m_explicit_abs_labels_shown;
    auto label = rdkit_extensions::get_atom_chirality_label(*m_atom, strip_abs);
    m_chirality_label_text = QString::fromStdString(label);

    auto chirality_label_position =
        findPositionInEmptySpace(true) * CHIRALITY_LABEL_DISTANCE_RATIO;
    m_chirality_label_rect =
        make_text_rect(m_fonts.m_chirality_fm, m_chirality_label_text);
    m_chirality_label_rect.moveCenter(chirality_label_position);
}

void AtomItem::clearLabels()
{
    for (auto string :
         {m_main_label_text, m_charge_and_radical_label_text,
          m_H_count_label_text, m_chirality_label_text, m_mapping_label_text}) {
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
        m_isotope_label_rect =
            make_text_rect(m_fonts.m_subscript_fm, m_isotope_label_text);
    }
}

void AtomItem::updateMappingLabel()
{
    auto mapping_n = m_atom->getAtomMapNum();
    // when no mapping is set RDKit returns 0
    if (mapping_n != 0) {
        m_mapping_label_text = QString("[%1]").arg(mapping_n);
        m_mapping_label_rect =
            make_text_rect(m_fonts.m_mapping_fm, m_mapping_label_text);
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
        m_charge_and_radical_label_rect = make_text_rect(
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
    m_isotope_label_rect.translate(
        QPointF(m_main_label_rect.left(), m_main_label_rect.center().y()) -
        m_isotope_label_rect.bottomRight());

    m_charge_and_radical_label_rect.translate(
        QPointF(m_main_label_rect.right(), m_main_label_rect.center().y()) -
        m_charge_and_radical_label_rect.bottomLeft());
    m_mapping_label_rect.translate(
        QPointF(m_main_label_rect.right(), m_main_label_rect.center().y()) -
        m_mapping_label_rect.bottomLeft());

    if (m_mapping_label_rect.isValid()) {
        m_charge_and_radical_label_rect.translate(m_mapping_label_rect.width(),
                                                  0);
    }
    label_rect = label_rect.united(m_isotope_label_rect);
    label_rect = label_rect.united(m_mapping_label_rect);

    auto H_direction = findHsDirection();
    // charge and radical text is next to the main label, on the right, unless
    // hydrogens are drawn to the right, in which case they come before charge
    // and radicals. (e.g. H3O+ but NH4+)
    if (H_direction != HsDirection::RIGHT) {
        label_rect = label_rect.united(m_charge_and_radical_label_rect);
    }
    const auto H_LABEL_VERTICAL_OFFSET =
        m_H_count_label_rect.height() * H_COUNT_LABEL_VERTICAL_OFFSET_RATIO;
    switch (H_direction) {
        case HsDirection::LEFT:
            m_H_count_label_rect.translate(
                label_rect.bottomLeft() -
                QPointF(m_H_count_label_rect.right(),
                        m_H_count_label_rect.top() + H_LABEL_VERTICAL_OFFSET));
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
                        m_H_count_label_rect.top() + H_LABEL_VERTICAL_OFFSET));
            label_rect.setRight(m_H_count_label_rect.right());
            m_charge_and_radical_label_rect.translate(
                m_H_count_label_rect.right() -
                    m_charge_and_radical_label_rect.left(),
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
                        m_H_count_label_rect.top() + H_LABEL_VERTICAL_OFFSET));
            QPointF offset(0, label_rect.top() - m_H_count_label_rect.bottom());

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
                        m_H_count_label_rect.top() + H_LABEL_VERTICAL_OFFSET));
            break;
        }
    }
}

bool AtomItem::determineLabelIsVisible() const
{
    // note that this method does not take m_valence_error_is_visible into
    // account.  If that value is true (and up-to-date), then the label should
    // be visible.
    if (m_settings.m_carbon_labels == CarbonLabels::ALL) {
        return true;
    }
    if (m_atom->getAtomicNum() != static_cast<unsigned int>(Element::C)) {
        return true;
    }
    if (m_atom->getIsotope() != 0) {
        return true;
    }
    if (m_atom->getAtomMapNum() != 0) {
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

QRectF AtomItem::getChiralityLabelRect() const
{
    return m_chirality_label_rect;
}
QString AtomItem::getChiralityLabelText() const
{
    return m_chirality_label_text;
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
    if (!m_squiggle_path.isEmpty()) {
        painter->strokePath(m_squiggle_path, m_squiggle_pen);
    }
    if (m_label_is_visible) {
        painter->save();
        painter->setFont(m_fonts.m_main_label_font);
        painter->setPen(m_pen);
        painter->drawText(m_main_label_rect, Qt::AlignCenter,
                          m_main_label_text);
        // usually drawText only paints when the rect is valid, but this appears
        // to not be the case for SVG (SKETCH-2233), so we need to check
        // explicitly
        if (m_H_label_rect.isValid()) {
            painter->drawText(m_H_label_rect, Qt::AlignCenter, "H");
        }
        painter->setFont(m_fonts.m_subscript_font);
        painter->drawText(m_H_count_label_rect, Qt::AlignCenter,
                          m_H_count_label_text);
        painter->drawText(m_charge_and_radical_label_rect, Qt::AlignCenter,
                          m_charge_and_radical_label_text);
        painter->drawText(m_isotope_label_rect, Qt::AlignCenter,
                          m_isotope_label_text);

        // atom mappings
        painter->setPen(m_pen);
        painter->setFont(m_fonts.m_mapping_font);
        painter->drawText(m_mapping_label_rect, Qt::AlignCenter,
                          m_mapping_label_text);
        painter->restore();
    }
    if (!m_query_label_text.isEmpty()) {
        painter->save();
        painter->setFont(m_fonts.m_element_list_font);
        painter->setPen(m_pen);
        painter->drawText(m_query_label_rect, Qt::AlignCenter,
                          m_query_label_text);
        auto element_list_line =
            QLineF(QPointF(0, 0), m_query_label_rect.center());
        trim_line_to_rect(element_list_line, m_query_label_rect);
        // trim the line on the atom's end. No label is drawn there, so we
        // define a small rectangle. A 4 pixel margin will be added by
        // trim_line_to_rect
        auto SMALL = 0.1;
        auto small_rect = QRectF(-SMALL, -SMALL, SMALL * 2, SMALL * 2);
        trim_line_to_rect(element_list_line, small_rect);
        painter->setPen(m_element_list_line_pen);
        painter->drawLine(element_list_line);
        painter->restore();
    }
    if (m_settings.m_stereo_labels_shown) {
        painter->save();
        painter->setPen(m_chirality_pen);
        painter->setFont(m_fonts.m_chirality_font);
        painter->drawText(m_chirality_label_rect, Qt::AlignCenter,
                          m_chirality_label_text);
        painter->restore();
    }
}

const RDKit::Atom* AtomItem::getAtom() const
{
    return m_atom;
}

QPointF AtomItem::findPositionInEmptySpace(bool avoid_subrects) const
{
    auto mol_positions = get_relative_positions_of_atom_neighbors(m_atom);
    std::vector<QPointF> qpositions;
    qpositions.reserve(mol_positions.size());
    for (auto pos : mol_positions) {
        qpositions.push_back(to_scene_xy(pos));
    }

    if (avoid_subrects) {
        for (auto rect : getSubrects()) {
            if (rect == m_main_label_rect)
                continue;
            qpositions.push_back(rect.topLeft());
            qpositions.push_back(rect.topRight());
            qpositions.push_back(rect.bottomLeft());
            qpositions.push_back(rect.bottomRight());
            qpositions.push_back(rect.center());
        }
    }
    return best_placing_around_origin(qpositions);
}

std::vector<QRectF> AtomItem::getLabelRects() const
{
    return {m_main_label_rect,      m_mapping_label_rect,
            m_isotope_label_rect,   m_charge_and_radical_label_rect,
            m_H_count_label_rect,   m_H_label_rect,
            m_chirality_label_rect, m_query_label_rect};
}

bool AtomItem::determineValenceErrorIsVisible() const
{
    return m_settings.m_valence_errors_shown && m_atom->hasValenceViolation();
}

} // namespace sketcher
} // namespace schrodinger
