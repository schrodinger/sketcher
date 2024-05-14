#pragma once

#include <rdkit/GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <QColor>

namespace schrodinger
{
namespace sketcher
{

/**
 * Which carbon atoms should be labeled (i.e. which carbon atoms get a "C"
 * drawn)
 */
enum class CarbonLabels {
    /// No labels for any carbons (unless there's some other reason to label
    /// that atom, e.g., it's an isolated atom or it has a non-zero charge)
    NONE,
    /// Only label carbons bound to exactly one non-hydrogen atom
    TERMINAL,
    /// Label all atoms
    ALL
};

/**
 * Coloring scheme corresponding to what is available in RDKit rendering
 */
enum class ColorScheme { DEFAULT, AVALON, CDK, DARK_MODE, BLACK_WHITE };

/**
 * An object used to store all settings relevant to the display of atoms
 */
class AtomDisplaySettings
{
  public:
    AtomDisplaySettings();

    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single AtomDisplaySettings object with all AtomItems.
    explicit AtomDisplaySettings(const AtomDisplaySettings&) = default;

    /**
     * Set a pre-defined RDKit color scheme
     * @param scheme color scheme to apply to atom items
     */
    void setColorScheme(ColorScheme scheme);

    /**
     * Set a color scheme where all atoms use the same color
     * @param color the color for all atoms
     */
    void setMonochromeColorScheme(QColor color);

    /**
     * @param atomic_number the atomic number of the atom
     * @return atom color based on the active color scheme
     */
    QColor getAtomColor(int atomic_number) const;

    /**
     * Which carbon atoms should be labeled
     */
    CarbonLabels m_carbon_labels = CarbonLabels::NONE;

    /**
     * Whether to display a yellow circle around atoms with valence errors
     */
    bool m_valence_errors_shown = true;

    /**
     * Whether to display stereochemistry
     * labels when present
     */
    bool m_stereo_labels_shown = true;

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

  private:
    RDKit::ColourPalette m_color_palette;
};
} // namespace sketcher
} // namespace schrodinger
