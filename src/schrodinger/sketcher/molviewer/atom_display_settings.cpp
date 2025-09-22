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
                                         QColor carbon_color)
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
        case ColorScheme::WHITE_BLACK:
            // RDKit doesn't actually have a White-Black color scheme, so we
            // make our own
            m_color_palette.clear();
            m_color_palette[-1] = get_rdkit_dark_mode_carbon_color();
    }
    if (carbon_color.isValid()) {
        m_color_palette[-1] =
            RDKit::DrawColour(carbon_color.redF(), carbon_color.greenF(),
                              carbon_color.blueF(), carbon_color.alphaF());
        if (m_color_palette.contains(static_cast<int>(Element::C))) {
            m_color_palette[static_cast<int>(Element::C)] = m_color_palette[-1];
        }
    }
    m_annotation_color =
        (scheme == ColorScheme::DARK_MODE || scheme == ColorScheme::WHITE_BLACK)
            ? ANNOTATION_COLOR_DARK
            : ANNOTATION_COLOR;
}

void AtomDisplaySettings::setMonochromeColorScheme(const QColor& color)
{
    // the color scheme's color is overridden by passing a value for
    // carbon_color, so the actual BLACK_WHITE color scheme will effectively be
    // ignored by this call (other than its lack of heteroatom colors)
    setColorScheme(ColorScheme::BLACK_WHITE, color);
}

QColor AtomDisplaySettings::getAtomColor(const int atomic_number) const
{
    auto color = m_color_palette.at(-1); // default
    auto it = m_color_palette.find(atomic_number);
    if (it != m_color_palette.end()) {
        color = it->second;
    }
    return QColor::fromRgbF(color.r, color.g, color.b, color.a);
}

void AtomDisplaySettings::setSquigglePenScale(const qreal scale)
{
    m_squiggle_pen_width = scale * BOND_DEFAULT_PEN_WIDTH;
}

RDKit::DrawColour get_rdkit_dark_mode_carbon_color()
{
    RDKit::ColourPalette palette;
    RDKit::assignDarkModePalette(palette);
    auto carbon = static_cast<int>(Element::C);
    carbon = palette.contains(carbon) ? carbon : -1;
    return palette[carbon];
}

} // namespace sketcher
} // namespace schrodinger
