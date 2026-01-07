#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"

#include <cmath>

#include <QPainter>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

namespace schrodinger
{
namespace sketcher
{

AbstractMonomerItem::AbstractMonomerItem(const RDKit::Atom* monomer,
                                         const Fonts& fonts,
                                         QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
}

void AbstractMonomerItem::paint(QPainter* painter,
                                const QStyleOptionGraphicsItem* option,
                                QWidget* widget)
{
    painter->save();
    painter->setPen(m_border_pen);
    painter->setBrush(m_border_brush);
    painter->drawPath(m_border_path);
    painter->setPen(m_main_label_pen);
    painter->setFont(m_main_label_font);
    painter->drawText(m_main_label_rect, Qt::AlignCenter, m_main_label_text);
    painter->restore();
}

void AbstractMonomerItem::setMonomerColors(const QColor& background_color,
                                           const QColor& outline_color,
                                           const QColor& font_color)
{
    m_border_brush.setColor(background_color);
    m_border_pen.setColor(outline_color);
    m_main_label_pen.setColor(font_color);
}

QString elide_text(const std::string& text)
{
    auto qtext = QString::fromStdString(text);
    if (text.length() > MAX_MONOMER_LABEL_LENGTH) {
        // shorten the string and add ellipses
        return qtext.left(MAX_MONOMER_LABEL_LENGTH - 1) +
               QString::fromUtf16(u"\u2026");
    } else {
        return qtext;
    }
}

std::string get_monomer_res_name(const RDKit::Atom* const monomer)
{
    bool is_smiles = false;
    monomer->getPropIfPresent(SMILES_MONOMER, is_smiles);
    if (is_smiles) {
        return SMILES_PLACEHOLDER_TEXT;
    }
    const auto* monomer_info = monomer->getMonomerInfo();
    const auto* res_info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info);
    if (res_info == nullptr) {
        return monomer_info->getName();
    }
    return res_info->getResidueName();
}

bool is_standard_nucleotide(const std::string& res_name)
{
    return NUCLEIC_ACID_COLOR_BY_RES_NAME.contains(res_name);
}

std::pair<qreal, qreal>
get_rect_size_to_fit_label(const QString& main_label_text,
                           const QFontMetricsF& fm, const qreal standard_width,
                           const qreal standard_height)
{
    auto tight_rect = fm.tightBoundingRect(main_label_text);
    auto min_width =
        tight_rect.width() + MONOMER_RECT_MINIMUM_HORIZONTAL_TEXT_PADDING;
    qreal border_width =
        std::max({standard_width, min_width, MONOMER_BORDER_MIN_SIZE});
    qreal border_height = std::max(
        {standard_height, tight_rect.height(), MONOMER_BORDER_MIN_SIZE});
    return {border_width, border_height};
}

/**
 * Determine the size of a diamond or ellipse required to fit the given label
 * text
 * @param main_label_text the text to size the shape for
 * @param fm a font metrics object for the font that will be used to display the
 * text
 * @param standard_width the minimum width of the shape
 * @param standard_height the minimum height of the shape
 * @param exponent the exponent of the shape's equation (1 for diamond, 2 for
 * ellipse).
 *
 * @note We determine the bounding rectangle of the text, and then solve the
 * equation of a diamond:
 *
 *   |x/x_radius| + |y/y_radius| = 1
 *
 * or an ellipse:
 *
 *   (x/x_radius)^2 + (y/y_radius)^2 = 1
 *
 * for one corner of that rectangle, which gives us the required size of the
 * outer shape (the diamond or ellipse). Initially, we fix the height of this
 * shape.  However, if this results in an aspect ratio less than
 * DIAMOND_AND_ELLIPSE_MIN_ASPECT_RATIO, then we increase the height slightly to
 * maintain the minimum aspect ratio.
 */
static std::pair<qreal, qreal> get_diamond_or_ellipse_size_to_fit_label(
    const QString& main_label_text, const QFontMetricsF& fm,
    const qreal standard_width, const qreal standard_height,
    const qreal exponent)
{
    auto tight_rect = fm.tightBoundingRect(main_label_text);
    auto text_width = tight_rect.width();
    auto text_height = tight_rect.height();
    auto min_height = standard_height;

    auto height_diff =
        std::pow(standard_height, exponent) - std::pow(text_height, exponent);
    auto clamped_height_diff = std::max(MONOMER_BORDER_MIN_SIZE, height_diff);
    auto min_width = text_width * standard_height /
                     std::pow(clamped_height_diff, 1.0 / exponent);
    if (standard_height < (DIAMOND_AND_ELLIPSE_MIN_ASPECT_RATIO * min_width)) {
        // the shape is too wide for the height, so we need to increase the
        // height
        auto min_width_to_the_n =
            std::pow(text_height / DIAMOND_AND_ELLIPSE_MIN_ASPECT_RATIO,
                     exponent) +
            std::pow(text_width, exponent);
        auto clamped_min_width_to_the_n =
            std::max(MONOMER_BORDER_MIN_SIZE, min_width_to_the_n);
        min_width = std::pow(clamped_min_width_to_the_n, 1.0 / exponent);
        min_height = DIAMOND_AND_ELLIPSE_MIN_ASPECT_RATIO * min_width;
    }

    qreal border_width =
        std::max({standard_width, min_width, MONOMER_BORDER_MIN_SIZE});
    qreal border_height =
        std::max({standard_height, min_height, MONOMER_BORDER_MIN_SIZE});
    return {border_width, border_height};
}

std::pair<qreal, qreal> get_diamond_size_to_fit_label(
    const QString& main_label_text, const QFontMetricsF& fm,
    const qreal standard_width, const qreal standard_height)
{
    return get_diamond_or_ellipse_size_to_fit_label(
        main_label_text, fm, standard_width, standard_height, 1.0);
}

std::pair<qreal, qreal> get_ellipse_size_to_fit_label(
    const QString& main_label_text, const QFontMetricsF& fm,
    const qreal standard_width, const qreal standard_height)
{
    return get_diamond_or_ellipse_size_to_fit_label(
        main_label_text, fm, standard_width, standard_height, 2.0);
}

std::tuple<qreal, QColor, const QFont, const QFontMetricsF*>
get_border_and_font_settings_for_nucleic_acid(const std::string& res_name,
                                              const Fonts& fonts,
                                              const bool use_base_font)
{
    if (is_standard_nucleotide(res_name)) {
        // for standard residues, we use the same font for both the base and the
        // backbone
        return {STANDARD_NA_BORDER_LINE_WIDTH, STANDARD_NA_BORDER_COLOR,
                fonts.m_main_label_font, &fonts.m_main_label_fm};
    } else if (use_base_font) {
        return {OTHER_NA_BORDER_LINE_WIDTH, OTHER_NA_BORDER_COLOR,
                fonts.m_other_nucleic_acid_base_font,
                &fonts.m_other_nucleic_acid_base_fm};
    } else {
        return {OTHER_NA_BORDER_LINE_WIDTH, OTHER_NA_BORDER_COLOR,
                fonts.m_other_nucleic_acid_sugar_and_phosphate_font,
                &fonts.m_other_nucleic_acid_sugar_and_phosphate_fm};
    }
}

void set_path_to_rect(QPainterPath& path, const QRectF& rect,
                      const qreal highlighting_thickness)
{
    auto expanded_rect =
        rect.adjusted(-highlighting_thickness, -highlighting_thickness,
                      highlighting_thickness, highlighting_thickness);
    path.clear();
    path.addRect(expanded_rect);
}

QRectF rect_expanded_by_half_pen_width(const QRectF& rect,
                                       const qreal pen_width)
{
    auto half_pen_width = pen_width / 2.0;
    return rect.adjusted(-half_pen_width, -half_pen_width, half_pen_width,
                         half_pen_width);
}

QColor get_color_for_monomer(
    const std::string& res_name, rdkit_extensions::ChainType chain_type,
    const std::unordered_map<std::string, QColor>& color_by_res_name,
    const QColor& default_color)
{
    auto color_find = color_by_res_name.find(res_name);
    if (color_find != color_by_res_name.end()) {
        return color_find->second;
    }

    // Check Monomer DB for natural analog
    auto& monomer_db = rdkit_extensions::MonomerDatabase::instance();
    auto natural_analog = monomer_db.getNaturalAnalog(res_name, chain_type);

    // If a valid natural analog exists, try to get its color
    if (natural_analog.has_value() && *natural_analog != res_name) {
        auto analog_color_find = color_by_res_name.find(*natural_analog);
        if (analog_color_find != color_by_res_name.end()) {
            return analog_color_find->second;
        }
    }

    return default_color;
}

} // namespace sketcher
} // namespace schrodinger
