#pragma once

#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
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
 * An object used to store all settings relevant to AtomItems
 */
class AtomItemSettings
{
  public:
    AtomItemSettings();

    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single AtomItemSettings object with all AtomItems.
    explicit AtomItemSettings(const AtomItemSettings&) = default;

    /**
     * @param scheme color scheme to apply to atom items
     */
    void setColorScheme(ColorScheme scheme);

    /**
     * @param atomic_number the atomic number of the atom
     * @return atom color based on the active color scheme
     */
    QColor getAtomColor(int atomic_number);

    /**
     * Which carbon atoms should be labeled
     */
    CarbonLabels m_carbon_labels = CarbonLabels::NONE;

    /**
     * Whether to display a yellow circle around atoms with valence errors
     */
    bool m_valence_errors_shown = true;

  private:
    RDKit::ColourPalette m_color_palette;
};
} // namespace sketcher
} // namespace schrodinger
