#pragma once

#include <rdkit/GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <QColor>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used to store all settings relevant to the display of atoms
 */
class SKETCHER_API AbstractAtomOrBondDisplaySettings
{
  public:
    AbstractAtomOrBondDisplaySettings() = default;
    explicit AbstractAtomOrBondDisplaySettings(
        const AbstractAtomOrBondDisplaySettings&) = default;
    virtual ~AbstractAtomOrBondDisplaySettings() = default;

    /**
     * Set a pre-defined RDKit color scheme
     * @param scheme color scheme to apply to items
     * @param carbon_color the color for carbon atoms.  If the monochrome color
     * scheme is being used, this color will apply to all atoms.  If the color
     * is invalid (the default), then the color scheme's default color will be
     * used (or black for monochrome color schemes).
     */
    virtual void setColorScheme(const ColorScheme& scheme,
                                const QColor& carbon_color = QColor());

    /**
     * Set a color scheme where all atoms use the same color
     * @param color the color for all atoms
     */
    void setMonochromeColorScheme(const QColor& color);

    /// The color of the annotation text
    QColor m_annotation_color;

  protected:
    /**
     * Update the given color palette so it reflects the specified color scheme
     * @param[in] scheme the color scheme
     * @param[in] carbon_color the color for carbon atoms.  If invalid, then the
     * color scheme's default color will be used.
     * @param[out] color_palette the color palette to update
     */
    void
    populatePaletteForColorScheme(const ColorScheme& scheme,
                                  const QColor& carbon_color,
                                  RDKit::ColourPalette& color_palette) const;

    /**
     * @return the color specified in the given color palette for an atom with
     * the specified atomic number
     */
    QColor
    getAtomColorFromPalette(const int atomic_number,
                            const RDKit::ColourPalette& color_palette) const;
};

/**
 * @return the light gray color that RDKit uses to color carbons when using the
 * dark mode palette
 */
RDKit::DrawColour get_rdkit_dark_mode_carbon_color();

} // namespace sketcher
} // namespace schrodinger
