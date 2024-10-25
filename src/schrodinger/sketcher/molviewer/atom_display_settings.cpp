#include "schrodinger/sketcher/molviewer/atom_display_settings.h"

namespace schrodinger
{
namespace sketcher
{

AtomDisplaySettings::AtomDisplaySettings()
{
    setColorScheme(ColorScheme::DEFAULT);
}

void AtomDisplaySettings::setColorScheme(ColorScheme scheme)
{
    switch (scheme) {
        case ColorScheme::DEFAULT:
            RDKit::assignDefaultPalette(m_color_palette);
            break;
        case ColorScheme::AVALON:
            RDKit::assignAvalonPalette(m_color_palette);
            break;
        case ColorScheme::CDK:
            RDKit::assignCDKPalette(m_color_palette);
            break;
        case ColorScheme::DARK_MODE:
            RDKit::assignDarkModePalette(m_color_palette);
            break;
        case ColorScheme::BLACK_WHITE:
            RDKit::assignBWPalette(m_color_palette);
            break;
    }
}

void AtomDisplaySettings::setMonochromeColorScheme(QColor color)
{
    m_color_palette.clear();
    m_color_palette[-1] =
        RDKit::DrawColour(color.redF(), color.greenF(), color.blueF());
}

QColor AtomDisplaySettings::getAtomColor(int atomic_number) const
{
    auto color = m_color_palette.at(-1); // default
    auto it = m_color_palette.find(atomic_number);
    if (it != m_color_palette.end()) {
        color = it->second;
    }
    return QColor::fromRgbF(color.r, color.g, color.b, color.a);
}

void AtomDisplaySettings::setSquigglePenScale(qreal scale)
{
    m_squiggle_pen_width = scale * BOND_DEFAULT_PEN_WIDTH;
}

} // namespace sketcher
} // namespace schrodinger
