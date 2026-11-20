#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/RDGeneral/types.h>

#include <QLineF>
#include <QPainter>
#include <QPainterPath>
#include <QPainterPathStroker>
#include <QPointF>
#include <QSizeF>
#include <QString>
#include <Qt>
#include <unordered_map>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"
#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/queries.h"

namespace schrodinger
{
namespace sketcher
{

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
    AbstractAtomOrMonomerItem(atom, parent),
    m_fonts(fonts),
    m_settings(settings)
{
    setZValue(static_cast<qreal>(ZOrder::ATOM));
    m_valence_error_pen = QPen(VALENCE_ERROR_BORDER_COLOR);
    m_valence_error_pen.setWidthF(VALENCE_ERROR_BORDER_WIDTH);
    m_valence_error_pen.setStyle(Qt::DotLine);
    m_valence_error_pen.setCapStyle(Qt::RoundCap);
    m_valence_error_brush = QBrush(settings.m_valence_error_area_color);

    // setup chirality pen
    m_chirality_pen = QPen(settings.m_annotation_color);

    // use the default atom color for the squiggle pen (used for attachment
    // points)
    m_squiggle_pen = QPen(settings.getAtomColor(-1));

    m_query_label_line_pen = QPen(QUERY_LABEL_LINE_COLOR);
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

void AtomItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();
    clearLabels();
    m_squiggle_pen.setWidthF(m_settings.m_squiggle_pen_width);

    bool needs_additional_labels, show_Hs;
    std::tie(m_main_label_text, m_squiggle_path, m_label_is_visible,
             m_valence_error_is_visible, needs_additional_labels, show_Hs,
             m_query_label_text) = determineLabelType();

    if (m_label_is_visible) {
        // if USER_COLOR is present override the color with USER_COLOR
        QColor user_color;
        m_atom->getPropIfPresent(USER_COLOR, user_color);
        m_pen.setColor(user_color.isValid()
                           ? user_color
                           : m_settings.getAtomColor(m_atom->getAtomicNum()));

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
            make_text_rect(m_fonts.m_query_label_fm, m_query_label_text);
        // find a position for the query label
        auto position =
            findPositionInEmptySpace(false) * QUERY_LABEL_DISTANCE_RATIO;
        // if the atom has no bonds the query label should be placed on
        // top (rather than to the right) of the atom so that it doesn't overlap
        // with the atom label
        if (m_atom->getDegree() == 0) {
            position = QPointF(0, -QUERY_LABEL_DISTANCE_RATIO * BOND_LENGTH *
                                      VIEW_SCALE);
        }
        // if the query label is a single character, or if it's above or below
        // the main label (as opposed to to the left or right), then we center
        // the label at this position
        m_query_label_rect.moveCenter(position);
        if (qAbs(position.y()) < qAbs(position.x()) &&
            m_query_label_text.size() > 1) {
            // If a multi-letter query label is off to the right or left of the
            // main label, then we want to left-align or right-align the query
            // label instead. That way, a long query label will always extend
            // away from the main label (or away from the bond if the main label
            // is hidden) instead of covering it.
            //
            // For the alignment, we want to (roughly) center the first or last
            // character of the string at the calculated position. That way, a
            // multi-character string will start at the same spot as a
            // single-character string (which will always be center-aligned).
            // For consistency in label positioning, we use a "standard" letter
            // width instead of grabbing the actual character from
            // m_query_label_extext. In other words, we just pick an arbitrary
            // letter and use that width since we don't want the label position
            // to vary with the width of the first/last character. This would be
            // particularly noticeable for SMARTS queries, as the first and last
            // characters of those are always double quotes, which are much
            // narrower than most letters. Using the double quote character
            // width would make it look like SMARTS queries were placed further
            // away from the main label than other types of queries.
            qreal char_width =
                m_fonts.m_query_label_fm.tightBoundingRect("A").width();
            qreal translate_x_by =
                (m_query_label_rect.width() - char_width) / 2;
            if (position.x() < 0) {
                translate_x_by *= -1;
            }
            m_query_label_rect.translate(translate_x_by, 0);
        }
    }
    if ((m_settings.m_stereo_labels_visibility != StereoLabels::NONE) &&
        m_settings.m_simplified_stereo_annotation.empty()) {
        updateChiralityLabel();
    }
    updateHighlightingPaths();
    m_shape = m_predictive_highlighting_path;
    m_bounding_rect = m_shape.boundingRect();
    // the annotation label must be positioned after all other labels and the
    // highlighting paths so that it can avoid them. Note that the annotation
    // label is intentionally omitted from the highlighting and halo paths.
    updateAnnotationLabel();
    if (m_annotation_label_rect.isValid()) {
        m_bounding_rect |= m_annotation_label_rect;
    }

    setToolTip(getTooltip());
}

void AtomItem::updateHighlightingPaths()
{
    m_selection_highlighting_path =
        get_selection_highlighting_path_for_atom(m_atom);
    m_predictive_highlighting_path =
        get_predictive_highlighting_path_for_atom(m_atom);
    const auto rects = getSubrects();
    if (rects.empty()) {
        // if this is an unlabeled carbon
        return;
    }
    const auto highlighting_paths = {
        std::make_pair(&m_predictive_highlighting_path, 0.0),
        std::make_pair(&m_selection_highlighting_path,
                       ATOM_HIGHLIGHTING_ADDITIONAL_SELECTION_RADIUS)};
    if (rects.size() == 1 && rects.front() == m_main_label_rect &&
        ((m_atom->getAtomicNum() > 0 &&
          m_main_label_text == QString::fromStdString(m_atom->getSymbol())) ||
         m_main_label_text.size() <= 1)) {
        // if the label is just the element name and nothing else, or if it's a
        // single character (e.g. "A" or "M" wildcards, "D" or "T" isotopes),
        // then we use a circular label
        const auto label_radius =
            QLineF(m_main_label_rect.center(), m_main_label_rect.bottomRight())
                .length() +
            TEXT_RECT_MARGIN;
        if (label_radius > ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS) {
            for (auto [path, padding] : highlighting_paths) {
                // if the default circle size is too small to cover the entire
                // text (e.g., a multi-letter element like "Cl"), enlarge the
                // circle just enough to cover
                path->clear();
                path->addEllipse(m_main_label_rect.center(),
                                 label_radius + padding,
                                 label_radius + padding);
            }
        }
        return;
    }
    // there are additional components to the label, so we use a rounded
    // rectangle large enough to cover everything
    QRectF label_rect;
    if (!m_label_is_visible) {
        // Offset labels must not exclude the atom's default highlighting area.
        label_rect = m_predictive_highlighting_path.boundingRect();
    }
    for (const auto& rect : rects) {
        label_rect |= rect;
    }
    label_rect.adjust(-ATOM_LABEL_MARGIN, -ATOM_LABEL_MARGIN, ATOM_LABEL_MARGIN,
                      ATOM_LABEL_MARGIN);
    // use the default rectangle rounding unless that's too big for this
    // rectangle, in which case use the largest rounding that can fit
    const auto corner_radius =
        std::min({ATOM_LABEL_HIGHLIGHTING_CORNER_RADIUS,
                  label_rect.width() * 0.5, label_rect.height() * 0.5});
    for (auto [path, padding] : highlighting_paths) {
        path->clear();
        path->addRoundedRect(
            label_rect.adjusted(-padding, -padding, padding, padding),
            corner_radius + padding, corner_radius + padding);
    }
}

QString AtomItem::advancedPropertiesSmarts() const
{
    auto props = read_properties_from_atom(m_atom);

    auto advanced_props = get_only_advanced_properties(
        std::dynamic_pointer_cast<AtomQueryProperties>(props));
    auto [advanced_properties_atom, maybe_enhanced_stereo] =
        create_atom_with_properties(advanced_props);
    auto smarts =
        QString::fromStdString(get_atom_smarts(advanced_properties_atom.get()));
    // remove square brackets and the #6& at the beginning
    smarts.chop(1);
    smarts.remove(0, 4);
    return smarts;
}

QString AtomItem::getQueryLabel() const
{
    QString query_text;
    auto props = read_properties_from_atom(m_atom);
    auto query_props = *static_cast<const AtomQueryProperties*>(props.get());

    // if this is an allowed or disallowed list with no added features, we can
    // use a simplified label, otherwise we need to show the SMARTS
    if ((query_props.query_type == QueryType::ALLOWED_LIST ||
         query_props.query_type == QueryType::NOT_ALLOWED_LIST) &&
        !query_props.hasPropertiesBeyondQueryType()) {
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
        query_text = m_fonts.m_query_label_fm.elidedText(
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
    bool has_user_set_label = false;

    // if there's a user-set display label, always use that.
    if (m_atom->hasProp(RDKit::common_properties::_displayLabel)) {
        main_label_text = m_atom->getProp<std::string>(
            RDKit::common_properties::_displayLabel);
        has_user_set_label = true;
    } else if (m_atom->getAtomicNum() !=
                   rdkit_extensions::DUMMY_ATOMIC_NUMBER &&
               m_atom->hasProp(RDKit::common_properties::atomLabel)) {

        // if there's no user-set, but an atomLabel is present, always
        // display that on non-dummy atoms. Query atoms are dealt with below
        main_label_text =
            m_atom->getProp<std::string>(RDKit::common_properties::atomLabel);
        has_user_set_label = true;

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
            has_user_set_label = true;
        } else if (is_dummy_atom_for_variable_attachment_bond(m_atom)) {

            // dummy atoms for variable attachment bonds aren't shown, but we
            // still need to set label_is_visible to false here so that the
            // bound BondItem knows not to make space for the "*" label
            label_is_visible = false;
        } else if (m_atom->hasQuery()) {
            // Wildcard or list query. This needs to be after R groups and
            // attachment points, since it would trigger for both but we need to
            // special-case them
            auto props = read_properties_from_atom(m_atom);
            auto query_props =
                *static_cast<const AtomQueryProperties*>(props.get());
            if (query_props.query_type == QueryType::ALLOWED_LIST ||
                query_props.query_type == QueryType::NOT_ALLOWED_LIST) {
                // for allowed and not allowed lists hide the main label
                label_is_visible = false;
                query_label_text = getQueryLabel();
            } else {
                main_label_text =
                    get_label_from_atom_query(query_props.wildcard)
                        .toStdString();
                needs_additional_labels = true;
                H_labels_are_visible = false;
            }
        } else {
            // Unrecognized dummy atom.  Display any user-set labels
            main_label_text = m_atom->getSymbol();
        }
    } else if (m_atom->hasQuery()) {
        label_is_visible = determineLabelIsVisible();
        main_label_text = m_atom->getSymbol();
        needs_additional_labels = true;
        // hide implicit Hs on queries
        H_labels_are_visible = false;
    } else {
        // a "normal" atom, i.e. an atom that represents an element
        valence_error_is_visible = determineValenceErrorIsVisible();
        label_is_visible =
            valence_error_is_visible || determineLabelIsVisible();
        if (label_is_visible) {
            main_label_text = m_atom->getSymbol();
            // if the label is set, mark deuterium and tritium with D and T
            if (m_settings.m_show_symbol_for_H_isotopes &&
                m_atom->getAtomicNum() == 1) {
                if (m_atom->getIsotope() == 2) {
                    main_label_text = "D";
                } else if (m_atom->getIsotope() == 3) {
                    main_label_text = "T";
                }
            }
            needs_additional_labels = true;
        }
    }
    if (m_atom->hasQuery() && !has_user_set_label) {
        // for smarts queries, display the smarts, unless we already have a
        // user set label (from the atomLabel or _displayLAbel property)
        auto props = read_properties_from_atom(m_atom);
        auto query_props =
            *static_cast<const AtomQueryProperties*>(props.get());
        if (query_props.query_type == QueryType::SMARTS) {
            query_label_text = getQueryLabel();
            label_is_visible = false;
        }
        if (query_label_text.isEmpty()) {
            // if advanced properties are present on the query, show them as
            // SMARTS
            query_label_text = advancedPropertiesSmarts();
        }
    }
    // if the atom has no bonds and we are showing the element list or smarts
    // label, show an asterisk as the main label
    if (!query_label_text.isEmpty() && m_atom->getDegree() == 0 &&
        !label_is_visible) {
        main_label_text = "*";
        label_is_visible = true;
    }

    // Append atom index to the main label when debug mode is enabled
    if (m_settings.m_show_atom_indices && label_is_visible) {
        main_label_text += ":" + std::to_string(m_atom->getIdx());
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
    m_chirality_label_rect = QRectF();
    auto strip_abs = !m_settings.m_explicit_abs_labels_shown;
    auto show_unknown_chirality =
        (m_settings.m_stereo_labels_visibility == StereoLabels::ALL);
    auto label =
        get_atom_chirality_label(*m_atom, strip_abs, show_unknown_chirality);
    m_chirality_label_text = QString::fromStdString(label);
    if (m_chirality_label_text.isEmpty()) {
        return;
    }
    // find a position for the chirality label before the rectangle is created,
    // otherwise it will be considered as a rectangle to avoid
    auto new_position = findPositionInEmptySpace(true);
    m_chirality_label_rect =
        make_text_rect(m_fonts.m_chirality_fm, m_chirality_label_text);
    // the center of the label is positioned r + k points away from the atom,
    // where r is radius of the circle around the label and k depends on the
    // bond length. This way the label remains at a fixed distance from the
    // atom, regardless of chirality font size
    float label_half_diagonal =
        std::sqrt(
            m_chirality_label_rect.width() * m_chirality_label_rect.width() +
            m_chirality_label_rect.height() * m_chirality_label_rect.height()) /
        2;
    float distance = BOND_LENGTH * VIEW_SCALE * CHIRALITY_LABEL_DISTANCE_RATIO +
                     label_half_diagonal;

    auto chirality_label_position =
        new_position / (BOND_LENGTH * VIEW_SCALE) * distance;
    m_chirality_label_rect.moveCenter(chirality_label_position);
}

/**
 * Find a position for an atom annotation label that doesn't overlap any of the
 * given obstacles. Candidate positions are generated by rotating the label
 * around the atom (starting from preferred_angle) and by moving it away from
 * the atom. The clear candidate that is closest to any of the anchor points is
 * returned, with rotation away from preferred_angle used as a tie-breaker.
 *
 * @param label_size The size of the label to position
 * @param preferred_angle The preferred direction for the label, in degrees,
 * using the QLineF::angle convention
 * @param max_distance The maximum distance between the atom and the nearest
 * edge of the label
 * @param anchor_points The points that the label should be as close as
 * possible to, in coordinates relative to the atom
 * @param obstacles Everything that the label must not overlap, in coordinates
 * relative to the atom
 * @return The label's bounding rectangle, in coordinates relative to the atom,
 * or std::nullopt if no clear position was found
 */
static std::optional<QRectF> find_annotation_label_rect(
    const QSizeF& label_size, const qreal preferred_angle,
    const qreal max_distance, const std::vector<QPointF>& anchor_points,
    const QPainterPath& obstacles)
{
    auto distance_to_anchors = [&anchor_points](const QRectF& rect) {
        qreal min_distance = std::numeric_limits<qreal>::max();
        for (const auto& anchor_point : anchor_points) {
            min_distance = std::min(min_distance,
                                    get_distance_to_rect(anchor_point, rect));
        }
        return min_distance;
    };
    struct Candidate {
        qreal cost;
        QRectF rect;
    };
    std::vector<Candidate> candidates;
    const int num_angles = 360 / ANNOTATION_LABEL_ANGLE_STEP;
    for (int i = 0; i < num_angles; ++i) {
        // alternate between rotating clockwise and counterclockwise, i.e. 0, 1,
        // -1, 2, -2, ..., so that the preferred direction is tried first
        const int num_steps = (i + 1) / 2 * (i % 2 ? 1 : -1);
        const auto direction =
            QLineF::fromPolar(1.0, preferred_angle +
                                       num_steps * ANNOTATION_LABEL_ANGLE_STEP)
                .p2();
        // how far the center of the label must be from the atom (in this
        // direction) for the near edge of the label to touch the atom
        const qreal center_offset =
            std::abs(direction.x()) * label_size.width() / 2 +
            std::abs(direction.y()) * label_size.height() / 2;
        const qreal angle_cost =
            std::abs(num_steps) * ANNOTATION_LABEL_ANGLE_STEP_PENALTY;
        for (qreal distance = 0.0; distance <= max_distance;
             distance += ANNOTATION_LABEL_DISTANCE_STEP) {
            QRectF rect(QPointF(), label_size);
            rect.moveCenter(direction * (distance + center_offset));
            candidates.push_back(
                {distance_to_anchors(rect) + angle_cost, rect});
        }
    }
    // use a stable sort so that ties are broken in favor of the preferred
    // direction
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [](const auto& a, const auto& b) { return a.cost < b.cost; });
    const auto padding = ANNOTATION_LABEL_PADDING;
    for (const auto& candidate : candidates) {
        if (!obstacles.intersects(candidate.rect.adjusted(-padding, -padding,
                                                          padding, padding))) {
            return candidate.rect;
        }
    }
    return std::nullopt;
}

void AtomItem::updateAnnotationLabel()
{
    m_annotation_label_text = QString();
    m_annotation_label_rect = QRectF();
    std::string annotation;
    if (!m_atom->getPropIfPresent(USER_ANNOTATION, annotation) ||
        annotation.empty()) {
        return;
    }
    m_annotation_label_text = QString::fromStdString(annotation);
    const auto label_size =
        make_text_rect(m_fonts.m_chirality_fm, m_annotation_label_text).size();
    const auto preferred_angle =
        QLineF(QPointF(0, 0), findPositionInEmptySpace(true)).angle();
    const qreal max_distance =
        BOND_LENGTH * VIEW_SCALE * ANNOTATION_LABEL_MAX_DISTANCE_RATIO;
    // the label's center can't be any further away than this
    const qreal search_radius =
        max_distance +
        QLineF(QPointF(0, 0), QPointF(label_size.width(), label_size.height()))
            .length() +
        ANNOTATION_LABEL_PADDING;
    const auto anchor_points = getAnnotationAnchorPoints();
    auto rect = find_annotation_label_rect(
        label_size, preferred_angle, max_distance, anchor_points,
        calcAnnotationObstacles(search_radius));
    if (!rect.has_value()) {
        // There's no empty space around the atom, so settle for not overlapping
        // this atom's other labels and highlighting. Note that the selection
        // highlighting path covers all other labels.
        rect = find_annotation_label_rect(label_size, preferred_angle,
                                          max_distance, anchor_points,
                                          m_selection_highlighting_path);
    }
    if (!rect.has_value()) {
        // If we still can't find a suitable location, just stick the label as
        // far away as possible in the preferred direction.  Note that this
        // should never actually happen since the highlighting path is much
        // smaller than the search area, but we include this fallback just in
        // case
        rect = QRectF(QPointF(), label_size);
        rect->moveCenter(QLineF::fromPolar(max_distance, preferred_angle).p2());
    }
    m_annotation_label_rect = *rect;
}

std::vector<QPointF> AtomItem::getAnnotationAnchorPoints() const
{
    std::vector<QPointF> anchor_points;
    for (const auto& rect :
         {m_main_label_rect, m_isotope_label_rect,
          m_charge_and_radical_label_rect, m_mapping_label_rect, m_H_label_rect,
          m_H_count_label_rect}) {
        if (rect.isValid()) {
            anchor_points.push_back(rect.center());
        }
    }
    if (anchor_points.empty()) {
        // the label isn't visible (e.g. an unlabeled carbon), so use the
        // atom's position
        anchor_points.emplace_back(0, 0);
    }
    return anchor_points;
}

QPainterPath AtomItem::calcAnnotationObstacles(const qreal search_radius) const
{
    // the selection highlighting path covers this atom's labels, so this
    // prevents the annotation from overlapping the labels or a halo
    QPainterPath obstacles = m_selection_highlighting_path;
    obstacles.setFillRule(Qt::WindingFill);

    const auto& mol = m_atom->getOwningMol();
    const auto& conf = mol.getConformer();
    const auto& origin = conf.getAtomPos(m_atom->getIdx());
    auto get_local_pos = [&conf, &origin](const RDKit::Atom* atom) {
        return to_scene_xy(conf.getAtomPos(atom->getIdx()) - origin);
    };
    const QRectF search_rect(-search_radius, -search_radius, 2 * search_radius,
                             2 * search_radius);

    // Avoid other atoms (and any halos around them). Note that we don't know
    // the size of the other atoms' labels, so we approximate them using the
    // default highlighting circle.
    for (const auto* atom : mol.atoms()) {
        if (atom == m_atom) {
            continue;
        }
        const auto pos = get_local_pos(atom);
        if (search_rect.contains(pos)) {
            obstacles.addEllipse(pos, ATOM_SELECTION_HIGHLIGHTING_RADIUS,
                                 ATOM_SELECTION_HIGHLIGHTING_RADIUS);
        }
    }

    // Avoid bonds. We use the width of the bond highlighting path so that the
    // label also avoids double bond lines, wedges, and bond halos.
    QPainterPath bond_lines;
    for (const auto* bond : mol.bonds()) {
        const QLineF line(get_local_pos(bond->getBeginAtom()),
                          get_local_pos(bond->getEndAtom()));
        // pad the line's bounding rect so that horizontal and vertical lines
        // don't produce an empty rect, which would never intersect anything
        const auto line_rect =
            QRectF(line.p1(), line.p2()).normalized().adjusted(-1, -1, 1, 1);
        if (search_rect.intersects(line_rect)) {
            bond_lines.moveTo(line.p1());
            bond_lines.lineTo(line.p2());
        }
    }
    QPainterPathStroker stroker;
    stroker.setWidth(2 * BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH);
    stroker.setCapStyle(Qt::FlatCap);
    obstacles.addPath(stroker.createStroke(bond_lines));
    return obstacles;
}

void AtomItem::clearLabels()
{
    for (auto* text : {&m_main_label_text, &m_charge_and_radical_label_text,
                       &m_H_count_label_text, &m_chirality_label_text,
                       &m_mapping_label_text, &m_isotope_label_text,
                       &m_query_label_text, &m_annotation_label_text}) {
        text->clear();
    }
    for (auto rect_member : getLabelRectMembers()) {
        this->*rect_member = QRectF();
    }
    m_annotation_label_rect = QRectF();
}

void AtomItem::updateIsotopeLabel()
{
    // clear the isotope label text and rect
    m_isotope_label_text = "";
    m_isotope_label_rect = QRectF();

    // when no isotope is set RDKit returns 0
    auto isotope = m_atom->getIsotope();
    // if showing D or T for deuterium and tritium, we're not showing the
    // isotope label
    if (m_settings.m_show_symbol_for_H_isotopes &&
        (m_atom->getAtomicNum() == 1 && (isotope == 2 || isotope == 3))) {
        return;
    }
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

QString AtomItem::getTooltip() const
{
    const QString STEREO_PREFIX = "Stereo: ";
    const QString QUERY_PREFIX = "Query: ";

    QStringList tooltip_parts;

    // Check for chirality
    if (!m_chirality_label_text.isEmpty()) {
        tooltip_parts.append(STEREO_PREFIX + m_chirality_label_text);
    }

    // Check if this is a query atom (wildcard or has query constraints)
    auto props = read_properties_from_atom(m_atom);

    if (auto* query_props =
            dynamic_cast<const AtomQueryProperties*>(props.get())) {
        // If there are additional constraints, show them
        if (!m_query_label_text.isEmpty()) {
            tooltip_parts.append(QUERY_PREFIX + m_query_label_text);
        } else if (query_props->query_type == QueryType::WILDCARD) {
            // For plain wildcards, provide descriptive tooltips
            static const std::unordered_map<AtomQuery, QString>
                wildcard_descriptions = {
                    {AtomQuery::A, "Any heavy atom"},
                    {AtomQuery::AH, "Any heavy atom or hydrogen"},
                    {AtomQuery::Q, "Any heteroatom"},
                    {AtomQuery::QH, "Any heteroatom or hydrogen"},
                    {AtomQuery::M, "Any metal"},
                    {AtomQuery::MH, "Any metal or hydrogen"},
                    {AtomQuery::X, "Any halogen"},
                    {AtomQuery::XH, "Any halogen or hydrogen"},
                };

            auto it = wildcard_descriptions.find(query_props->wildcard);
            if (it != wildcard_descriptions.end()) {
                tooltip_parts.append(QUERY_PREFIX + it->second);
            }
        }
    }

    return tooltip_parts.join("\n");
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
    if (m_settings.m_show_atom_indices) {
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

const std::vector<QRectF> AtomItem::getSubrects() const
{
    std::vector<QRectF> subrects;
    for (auto rect : getLabelRects()) {
        if (rect.isValid()) {
            subrects.push_back(rect);
        }
    }
    return subrects;
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
        painter->setFont(m_fonts.m_query_label_font);
        painter->setPen(m_pen);
        painter->drawText(m_query_label_rect, Qt::AlignCenter,
                          m_query_label_text);
        if (!m_label_is_visible) {
            // if the main label is hidden (which happens for allowed/disallowed
            // lists and SMARTS queries), draw a small line between the query
            // label and the atom's center
            auto query_label_line =
                QLineF(QPointF(0, 0), m_query_label_rect.center());
            trim_line_to_rect(query_label_line, m_query_label_rect);
            // trim the line on the atom's end by defining a small rectangle.
            // A 4 pixel margin will be added by trim_line_to_rect
            auto SMALL = 0.1;
            auto small_rect = QRectF(-SMALL, -SMALL, SMALL * 2, SMALL * 2);
            trim_line_to_rect(query_label_line, small_rect);
            painter->setPen(m_query_label_line_pen);
            painter->drawLine(query_label_line);
        }
        painter->restore();
    }
    if (m_settings.m_stereo_labels_visibility != StereoLabels::NONE) {
        painter->save();
        painter->setPen(m_chirality_pen);
        painter->setFont(m_fonts.m_chirality_font);
        painter->drawText(m_chirality_label_rect, Qt::AlignCenter,
                          m_chirality_label_text);
        painter->restore();
    }
    // drawText usually only paints when the rect is valid, but this isn't the
    // case for SVG (SKETCH-2233), so we check explicitly
    if (m_annotation_label_rect.isValid()) {
        painter->save();
        painter->setPen(m_chirality_pen);
        painter->setFont(m_fonts.m_chirality_font);
        painter->drawText(m_annotation_label_rect, Qt::AlignCenter,
                          m_annotation_label_text);
        painter->restore();
    }
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

const std::vector<QRectF AtomItem::*>& AtomItem::getLabelRectMembers()
{
    static const std::vector<QRectF AtomItem::*> LABEL_RECT_MEMBERS = {
        &AtomItem::m_main_label_rect,
        &AtomItem::m_mapping_label_rect,
        &AtomItem::m_isotope_label_rect,
        &AtomItem::m_charge_and_radical_label_rect,
        &AtomItem::m_H_count_label_rect,
        &AtomItem::m_H_label_rect,
        &AtomItem::m_chirality_label_rect,
        &AtomItem::m_query_label_rect};
    return LABEL_RECT_MEMBERS;
}

std::vector<QRectF> AtomItem::getLabelRects() const
{
    std::vector<QRectF> rects;
    for (auto rect_member : getLabelRectMembers()) {
        rects.push_back(this->*rect_member);
    }
    return rects;
}

bool AtomItem::determineValenceErrorIsVisible() const
{
    return m_settings.m_valence_errors_shown && m_atom->hasValenceViolation();
}

} // namespace sketcher
} // namespace schrodinger
