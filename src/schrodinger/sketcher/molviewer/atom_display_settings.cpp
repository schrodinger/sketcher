#include "schrodinger/sketcher/molviewer/atom_display_settings.h"

#include <cstdlib>

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

AtomDisplaySettings::AtomDisplaySettings()
{
    setColorScheme(ColorScheme::DEFAULT);

    // Check if debug mode should be enabled via environment variable
    const char* debug_mode_env = std::getenv("SCHRODINGER_SKETCHER_DEBUG_MODE");
    if (debug_mode_env != nullptr && debug_mode_env[0] != '\0') {
        // Any non-empty value enables the feature
        m_show_atom_indices = true;
    }
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
