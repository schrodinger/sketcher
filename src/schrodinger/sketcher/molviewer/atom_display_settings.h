#pragma once

#include <rdkit/GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <QColor>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used to store all settings relevant to the display of atoms
 */
class SKETCHER_API AtomDisplaySettings
    : public AbstractAtomOrBondDisplaySettings
{
  public:
    AtomDisplaySettings();

    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single AtomDisplaySettings object with all AtomItems.
    explicit AtomDisplaySettings(const AtomDisplaySettings&) = default;

    virtual ~AtomDisplaySettings() = default;

    /**
     * Set a pre-defined RDKit color scheme
     * @param scheme color scheme to apply to atom items
     * @param carbon_color the color for carbon atoms.  If the monochrome color
     * scheme is being used, this color will apply to all atoms.  If the color
     * is invalid (the default), then the color scheme's default color will be
     * used (or black for monochrome color schemes).
     */
    void setColorScheme(const ColorScheme& scheme,
                        const QColor& carbon_color = QColor()) override;

    /**
     * @param atomic_number the atomic number of the atom
     * @return atom color based on the active color scheme
     */
    QColor getAtomColor(const int atomic_number) const;

    /**
     * Scale the width of the pen used to draw attachment point squiggles
     * @param scale The scale to use.  A value of 1.0 will result in the default
     * pen width.
     */
    void setSquigglePenScale(const qreal scale);

    /**
     * Which carbon atoms should be labeled
     */
    CarbonLabels m_carbon_labels = CarbonLabels::NONE;

    /**
     * Whether to display a yellow circle around atoms with valence errors
     */
    bool m_valence_errors_shown = true;

    /**
     * Whether to display stereochemistry labels when present
     */
    StereoLabels m_stereo_labels_visibility = StereoLabels::ALL;

    /**
     * Whether to display the word "abs" in enhanced stereo absolute
     * stereochemistry labels
     */
    bool m_explicit_abs_labels_shown = false;

    /**
     * Whether to display a simplified stereochemistry annotation for a compound
     * with only one ehnaced stereo group that is either "or" or "and". If this
     * is true a label like "or enantiomer" or "end enantiomer" will be
     * displayed insted of the atom labels.
     */
    bool m_show_simplified_stereo_annotation = false;

    /**
     * whether to display H isotopes with a symbol (D for deuterium and T for
     * tritium)
     */
    bool m_show_symbol_for_H_isotopes = false;

    /**
     * Whether to display RDKit mol atom indices as labels (developer option)
     */
    bool m_show_atom_indices = false;

    /**
     * The width of the pen used to draw attachment point squiggles.
     */
    qreal m_squiggle_pen_width = BOND_DEFAULT_PEN_WIDTH;

    QColor m_valence_error_area_color = VALENCE_ERROR_AREA_COLOR;

    /**
     * The simplified stereo annotation to display, if any.  This is set by the
     * Scene when the display settings are changed and it is used by AtomItems
     * to decide whether to show or hide their atomic labels.
     */
    std::string m_simplified_stereo_annotation = "";

  private:
    RDKit::ColourPalette m_color_palette;
};

} // namespace sketcher
} // namespace schrodinger
