/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

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
 * An object used to store all settings relevant to AtomItems
 */
class AtomItemSettings
{
  public:
    AtomItemSettings() = default;
    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single AtomItemSettings object with all AtomItems.
    explicit AtomItemSettings(const AtomItemSettings&) = default;
    /// Which carbon atoms should be labeled
    CarbonLabels m_carbon_labels = CarbonLabels::NONE;
    /// Whether to display a yellow circle around atoms with valence errors
    bool m_valence_errors_shown = true;
};
} // namespace sketcher
} // namespace schrodinger
