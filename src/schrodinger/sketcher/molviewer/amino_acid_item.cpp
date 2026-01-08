#include "schrodinger/sketcher/molviewer/amino_acid_item.h"

#include <cctype>

#include <QPainter>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

namespace schrodinger
{
namespace sketcher
{

enum class AminoAcidType { STANDARD, D, OTHER };

AminoAcidItem::AminoAcidItem(const RDKit::Atom* monomer, const Fonts& fonts,
                             QGraphicsItem* parent) :
    AbstractMonomerItem(monomer, fonts, parent)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int AminoAcidItem::type() const
{
    return Type;
}

static AminoAcidType get_amino_acid_type(const std::string& res_name)
{
    if (AMINO_ACID_COLOR_BY_RES_NAME.contains(res_name)) {
        return AminoAcidType::STANDARD;
    } else if (std::tolower(res_name[0]) == 'd') {
        return AminoAcidType::D;
    } else {
        return AminoAcidType::OTHER;
    }
}

/**
 * Update the given path so it contains only a rounded rectangle of the
 * specified size, using the amino acid rounding radius
 * @param[in,out] path the path to modify
 * @param[in] rect the rectangle to add to the path
 * @param[in] highlighting_thickness additional size to add to the rectangle.
 * This option is intended for generating selection and predictive highlighting
 * paths.
 */
static void set_path_to_rounded_rect(QPainterPath& highlighting_path,
                                     const QRectF& border_rect,
                                     const qreal highlighting_thickness = 0)
{
    auto selection_highlighting_rect =
        border_rect.adjusted(-highlighting_thickness, -highlighting_thickness,
                             highlighting_thickness, highlighting_thickness);
    highlighting_path.clear();
    highlighting_path.addRoundedRect(
        selection_highlighting_rect,
        AA_ROUNDING_RADIUS + highlighting_thickness,
        AA_ROUNDING_RADIUS + highlighting_thickness);
}

void AminoAcidItem::updateCachedData()
{
    prepareGeometryChange();

    auto res_name = get_monomer_res_name(m_atom);
    m_main_label_text = elide_text(res_name);
    qreal standard_border_width, standard_border_height, border_line_width;
    QColor border_color;
    const QFontMetricsF* main_label_fm = nullptr;
    switch (get_amino_acid_type(res_name)) {
        case AminoAcidType::STANDARD:
            standard_border_width = STANDARD_AA_BORDER_WIDTH;
            standard_border_height = STANDARD_AA_BORDER_HEIGHT;
            border_line_width = STANDARD_AA_BORDER_LINE_WIDTH;
            border_color = STANDARD_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_main_label_font;
            main_label_fm = &m_fonts.m_main_label_fm;
            break;
        case AminoAcidType::D:
            standard_border_width = D_AA_BORDER_WIDTH;
            standard_border_height = D_AA_BORDER_HEIGHT;
            border_line_width = D_AA_BORDER_LINE_WIDTH;
            border_color = D_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_d_amino_acid_font;
            main_label_fm = &m_fonts.m_d_amino_acid_fm;
            break;
        case AminoAcidType::OTHER:
            standard_border_width = OTHER_AA_BORDER_WIDTH;
            standard_border_height = OTHER_AA_BORDER_HEIGHT;
            border_line_width = OTHER_AA_BORDER_LINE_WIDTH;
            border_color = OTHER_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_other_amino_acid_font;
            main_label_fm = &m_fonts.m_other_amino_acid_fm;
            break;
    }
    auto [border_width, border_height] = get_rect_size_to_fit_label(
        m_main_label_text, *main_label_fm, standard_border_width,
        standard_border_height);
    m_border_pen.setWidthF(border_line_width);
    m_border_pen.setColor(border_color);
    auto border_rect = QRectF(-border_width / 2, -border_height / 2,
                              border_width, border_height);
    set_path_to_rounded_rect(m_border_path, border_rect);
    set_path_to_rounded_rect(m_selection_highlighting_path, border_rect,
                             MONOMER_SELECTION_HIGHLIGHTING_THICKNESS);
    set_path_to_rounded_rect(m_predictive_highlighting_path, border_rect,
                             MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS);
    m_main_label_rect = border_rect;
    m_bounding_rect =
        rect_expanded_by_half_pen_width(border_rect, border_line_width);
    set_path_to_rounded_rect(m_shape, border_rect, border_line_width / 2.0);
    auto rect_color = get_color_for_monomer(
        res_name, rdkit_extensions::ChainType::PEPTIDE,
        AMINO_ACID_COLOR_BY_RES_NAME, DEFAULT_AA_BACKGROUND_COLOR);
    m_border_brush.setColor(rect_color);
}

} // namespace sketcher
} // namespace schrodinger
