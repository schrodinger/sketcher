#include "schrodinger/sketcher/molviewer/atom_item_settings.h"

namespace schrodinger
{
namespace sketcher
{

AtomItemSettings::AtomItemSettings()
{
    setColorScheme(ColorScheme::DEFAULT);
}

void AtomItemSettings::setColorScheme(ColorScheme scheme)
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

QColor AtomItemSettings::getAtomColor(int atomic_number)
{
    auto color = m_color_palette[-1]; // default
    auto it = m_color_palette.find(atomic_number);
    if (it != m_color_palette.end()) {
        color = it->second;
    }
    return QColor::fromRgbF(color.r, color.g, color.b, color.a);
}

} // namespace sketcher
} // namespace schrodinger
