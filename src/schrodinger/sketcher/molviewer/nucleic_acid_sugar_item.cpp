#include "schrodinger/sketcher/molviewer/nucleic_acid_sugar_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidSugarItem::NucleicAcidSugarItem(const RDKit::Atom* monomer,
                                           const Fonts& fonts,
                                           QGraphicsItem* parent) :
    AbstractMonomerItem(monomer, fonts, parent)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    m_border_brush.setColor(NA_BACKBONE_COLOR);
    updateCachedData();
}

int NucleicAcidSugarItem::type() const
{
    return Type;
}

void NucleicAcidSugarItem::updateCachedData()
{
    prepareGeometryChange();

    auto res_name = get_monomer_res_name(m_atom);
    m_main_label_text = elide_text(res_name);

    auto [border_pen_width, border_color, border_color_dark_bg, main_label_font,
          main_label_fm] =
        get_border_and_font_settings_for_nucleic_acid(res_name, m_fonts);
    m_main_label_font = main_label_font;
    m_border_pen.setWidthF(border_pen_width);
    m_border_color = border_color;
    m_border_color_dark_bg = border_color_dark_bg;
    m_border_pen.setColor(getBorderColor());
    auto [border_width, border_height] = get_rect_size_to_fit_label(
        m_main_label_text, *main_label_fm,
        scaleBasedOnFontSize(NA_SUGAR_BORDER_WIDTH),
        scaleBasedOnFontSize(NA_SUGAR_BORDER_HEIGHT));

    auto border_rect = QRectF(-border_width / 2, -border_height / 2,
                              border_width, border_height);
    set_path_to_rect(m_border_path, border_rect);
    set_path_to_rect(m_selection_highlighting_path, border_rect,
                     MONOMER_SELECTION_HIGHLIGHTING_THICKNESS);
    set_path_to_rect(m_predictive_highlighting_path, border_rect,
                     MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS);
    m_main_label_left_baseline =
        center_text_in_rect(m_main_label_text, *main_label_fm, border_rect);
    m_bounding_rect =
        rect_expanded_by_half_pen_width(border_rect, border_pen_width);
    set_path_to_rect(m_shape, border_rect, border_pen_width / 2.0);
}

} // namespace sketcher
} // namespace schrodinger
