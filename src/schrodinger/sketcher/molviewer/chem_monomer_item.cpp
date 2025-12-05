#include "schrodinger/sketcher/molviewer/chem_monomer_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

ChemMonomerItem::ChemMonomerItem(const RDKit::Atom* monomer, const Fonts& fonts,
                                 QGraphicsItem* parent) :
    AbstractMonomerItem(monomer, fonts, parent)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    m_border_brush.setColor(CHEM_MONOMER_RECT_COLOR);
    m_main_label_font = fonts.m_main_label_font;
    m_border_pen.setWidth(CHEM_MONOMER_BORDER_LINE_WIDTH);
    m_border_pen.setColor(CHEM_MONOMER_BORDER_COLOR);
    updateCachedData();
}

int ChemMonomerItem::type() const
{
    return Type;
}

void ChemMonomerItem::updateCachedData()
{
    prepareGeometryChange();
    auto res_name = get_monomer_res_name(m_atom);
    m_main_label_text = elide_text(res_name);
    auto [border_width, border_height] = get_rect_size_to_fit_label(
        m_main_label_text, m_fonts.m_main_label_fm, CHEM_MONOMER_BORDER_WIDTH,
        CHEM_MONOMER_BORDER_HEIGHT);
    auto border_rect = QRectF(-border_width / 2, -border_height / 2,
                              border_width, border_height);
    set_path_to_rect(m_border_path, border_rect);
    set_path_to_rect(m_selection_highlighting_path, border_rect,
                     MONOMER_SELECTION_HIGHLIGHTING_THICKNESS);
    set_path_to_rect(m_predictive_highlighting_path, border_rect,
                     MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS);
    m_main_label_rect = border_rect;
    // the bounding rect needs to account for the thickness of the border
    // outline (which is centered on the outline, so half of it is *outside*
    // border_rect)
    auto border_line_half_width = CHEM_MONOMER_BORDER_LINE_WIDTH / 2.0;
    m_bounding_rect =
        border_rect.adjusted(-border_line_half_width, -border_line_half_width,
                             border_line_half_width, border_line_half_width);
    m_shape = m_border_path;
}

} // namespace sketcher
} // namespace schrodinger
