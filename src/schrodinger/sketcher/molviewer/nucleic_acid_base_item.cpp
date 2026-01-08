#include "schrodinger/sketcher/molviewer/nucleic_acid_base_item.h"

#include <cmath>

#include <QPainter>
#include <QPointF>
#include <QPolygonF>

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

NucleicAcidBaseItem::NucleicAcidBaseItem(const RDKit::Atom* monomer,
                                         const Fonts& fonts,
                                         QGraphicsItem* parent) :
    AbstractMonomerItem(monomer, fonts, parent)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));

    updateCachedData();
}

int NucleicAcidBaseItem::type() const
{
    return Type;
}

/**
 * Update the given path so it contains only a diamond of the specified size.
 * @param[in,out] path the path to modify
 * @param[in] width the width of the diamond to create
 * @param[in] height the height of the diamond to create
 * @param[in] highlighting_thickness additional height to add to the diamond for
 * purposes of selection and predictive highlighting. Additional width will be
 * added as well to maintain the aspect ratio of the diamond.
 */
static void set_path_to_diamond(QPainterPath& path, const qreal width,
                                const qreal height,
                                qreal highlighting_thickness = 0.0)
{
    // the passed-in highlighting_thickness should specify the thickness as
    // measured perpendicular to the diamond (so that the size of the diamond
    // highlighting matches the size of the square highlighting), so calculate
    // the corresponding axis-aligned thickness
    highlighting_thickness *= std::sqrt(2.0);
    // scale the additional distance so that the shape maintains the same aspect
    // ratio as the original
    auto half_width = width / 2 + highlighting_thickness * width / height;
    auto half_height = height / 2 + highlighting_thickness;
    path.clear();
    auto polygon = QPolygonF({
        QPointF(-half_width, 0),
        QPointF(0, half_height),
        QPointF(half_width, 0),
        QPointF(0, -half_height),
    });
    path.addPolygon(polygon);
    path.closeSubpath();
}

void NucleicAcidBaseItem::updateCachedData()
{
    prepareGeometryChange();

    auto res_name = get_monomer_res_name(m_atom);
    m_main_label_text = elide_text(res_name);

    auto [border_pen_width, border_color, main_label_font, main_label_fm] =
        get_border_and_font_settings_for_nucleic_acid(
            res_name, m_fonts,
            /* use_base_font = */ true);
    m_main_label_font = main_label_font;
    m_border_pen.setWidthF(border_pen_width);
    m_border_pen.setColor(border_color);

    auto [border_width, border_height] = get_diamond_size_to_fit_label(
        m_main_label_text, *main_label_fm, NA_BASE_BORDER_WIDTH,
        NA_BASE_BORDER_HEIGHT);
    set_path_to_diamond(m_border_path, border_width, border_height);
    set_path_to_diamond(m_selection_highlighting_path, border_width,
                        border_height,
                        MONOMER_SELECTION_HIGHLIGHTING_THICKNESS);
    set_path_to_diamond(m_predictive_highlighting_path, border_width,
                        border_height,
                        MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS);
    m_main_label_rect = m_border_path.boundingRect();
    set_path_to_diamond(m_shape, border_width, border_height,
                        border_pen_width / 2.0);
    m_bounding_rect = m_shape.boundingRect();

    // According to HELM, DNA is a subtype of RNA, so both use ChainType::RNA
    auto rect_color = get_color_for_monomer(
        res_name, rdkit_extensions::ChainType::RNA,
        NUCLEIC_ACID_COLOR_BY_RES_NAME, DEFAULT_NA_BACKGROUND_COLOR);
    m_border_brush.setColor(rect_color);
}

} // namespace sketcher
} // namespace schrodinger
