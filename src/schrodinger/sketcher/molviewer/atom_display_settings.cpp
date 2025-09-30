#include "schrodinger/sketcher/molviewer/atom_display_settings.h"

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

AtomDisplaySettings::AtomDisplaySettings()
{
    setColorScheme(ColorScheme::DEFAULT);
}

void AtomDisplaySettings::setColorScheme(const ColorScheme& scheme,
                                         const QColor& carbon_color)
{
    // update the annotation color
    AbstractAtomOrBondDisplaySettings::setColorScheme(scheme, carbon_color);
    m_valence_error_area_color =
        (scheme == ColorScheme::DARK_MODE || scheme == ColorScheme::WHITE_BLACK)
            ? VALENCE_ERROR_AREA_COLOR_DARK_BG
            : VALENCE_ERROR_AREA_COLOR;
    populatePaletteForColorScheme(scheme, carbon_color, m_color_palette);
}

QColor AtomDisplaySettings::getAtomColor(const int atomic_number) const
{
    return getAtomColorFromPalette(atomic_number, m_color_palette);
}

void AtomDisplaySettings::setSquigglePenScale(const qreal scale)
{
    m_squiggle_pen_width = scale * BOND_DEFAULT_PEN_WIDTH;
}

} // namespace sketcher
} // namespace schrodinger
