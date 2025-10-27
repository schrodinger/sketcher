#include "schrodinger/sketcher/molviewer/abstract_atom_or_bond_display_settings.h"

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

void AbstractAtomOrBondDisplaySettings::setColorScheme(
    const ColorScheme& scheme, const QColor& carbon_color)
{
    m_annotation_color = DARK_MODE_COLOR_SCHEMES.contains(scheme)
                             ? ANNOTATION_COLOR_DARK_BG
                             : ANNOTATION_COLOR;
}

void AbstractAtomOrBondDisplaySettings::populatePaletteForColorScheme(
    const ColorScheme& scheme, const QColor& carbon_color,
    RDKit::ColourPalette& color_palette) const
{
    switch (scheme) {
        case ColorScheme::DEFAULT:
            RDKit::assignDefaultPalette(color_palette);
            break;
        case ColorScheme::AVALON:
            RDKit::assignAvalonPalette(color_palette);
            break;
        case ColorScheme::CDK:
            RDKit::assignCDKPalette(color_palette);
            break;
        case ColorScheme::DARK_MODE:
            RDKit::assignDarkModePalette(color_palette);
            break;
        case ColorScheme::BLACK_WHITE:
            RDKit::assignBWPalette(color_palette);
            break;
        case ColorScheme::WHITE_BLACK:
            // RDKit doesn't actually have a White-Black color scheme, so we
            // make our own
            color_palette.clear();
            color_palette[-1] = get_rdkit_dark_mode_carbon_color();
    }
    if (carbon_color.isValid()) {
        color_palette[-1] =
            RDKit::DrawColour(carbon_color.redF(), carbon_color.greenF(),
                              carbon_color.blueF(), carbon_color.alphaF());
        if (color_palette.contains(static_cast<int>(Element::C))) {
            color_palette[static_cast<int>(Element::C)] = color_palette[-1];
        }
    }
}

void AbstractAtomOrBondDisplaySettings::setMonochromeColorScheme(
    const QColor& color)
{
    // the color scheme's color is overridden by passing a value for
    // carbon_color, so the actual BLACK_WHITE color scheme will effectively be
    // ignored by this call (other than its lack of heteroatom colors)
    setColorScheme(ColorScheme::BLACK_WHITE, color);
}

QColor AbstractAtomOrBondDisplaySettings::getAtomColorFromPalette(
    const int atomic_number, const RDKit::ColourPalette& color_palette) const
{
    auto color = color_palette.at(-1); // default
    auto it = color_palette.find(atomic_number);
    if (it != color_palette.end()) {
        color = it->second;
    }
    return QColor::fromRgbF(color.r, color.g, color.b, color.a);
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
