#include "schrodinger/sketcher/molviewer/nucleic_acid_phosphate_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidPhosphateItem::NucleicAcidPhosphateItem(const RDKit::Atom* monomer,
                                                   const Fonts& fonts,
                                                   QGraphicsItem* parent) :
    AbstractMonomerItem(monomer, fonts, parent)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    m_border_brush.setColor(NA_BACKBONE_COLOR);
    updateCachedData();
}

int NucleicAcidPhosphateItem::type() const
{
    return Type;
}

/**
 * Update the given path so it contains only an ellipse of the specified size
 * @param[in,out] path the path to modify
 * @param[in] rect the bounding rectangle of the ellipse
 * @param[in] highlighting_thickness additional distance to add to the radii of
 * the ellipse. This option is intended for generating selection and predictive
 * highlighting paths.
 */
static void set_path_to_ellipse(QPainterPath& path, const QRectF& rect,
                                const qreal highlighting_thickness = 0)
{
    auto expanded_rect =
        rect.adjusted(-highlighting_thickness, -highlighting_thickness,
                      highlighting_thickness, highlighting_thickness);
    path.clear();
    path.addEllipse(expanded_rect);
}

void NucleicAcidPhosphateItem::updateCachedData()
{
    prepareGeometryChange();

    auto res_name = get_monomer_res_name(m_atom);
    m_main_label_text = elide_text(res_name);

    auto [border_pen_width, border_color, main_label_font, main_label_fm] =
        get_border_and_font_settings_for_nucleic_acid(res_name, m_fonts);
    m_main_label_font = main_label_font;
    m_border_pen.setWidthF(border_pen_width);
    m_border_pen.setColor(border_color);

    auto [border_width, border_height] = get_ellipse_size_to_fit_label(
        m_main_label_text, *main_label_fm, NA_PHOSPHATE_BORDER_WIDTH,
        NA_PHOSPHATE_BORDER_HEIGHT);

    auto border_rect = QRectF(-border_width / 2, -border_height / 2,
                              border_width, border_height);
    set_path_to_ellipse(m_border_path, border_rect);
    set_path_to_ellipse(m_selection_highlighting_path, border_rect,
                        MONOMER_SELECTION_HIGHLIGHTING_THICKNESS);
    set_path_to_ellipse(m_predictive_highlighting_path, border_rect,
                        MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS);
    m_main_label_rect = border_rect;
    m_bounding_rect =
        rect_expanded_by_half_pen_width(border_rect, border_pen_width);
    set_path_to_rect(m_shape, border_rect, border_pen_width / 2.0);
}

} // namespace sketcher
} // namespace schrodinger
