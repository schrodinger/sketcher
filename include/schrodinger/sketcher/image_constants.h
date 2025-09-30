#pragma once

#include <unordered_set>

#include <QColor>

namespace schrodinger
{
namespace sketcher
{

/**
 * Enumerate the conditions under which heteroatoms should be colored.
 *
 * This enum should correspond directly to the Maestro preference enum
 * MMENUM_LOOKUP_2D_COLOR.
 */
enum class ColorHeteroatomsMode { NEVER, LIGHT_BACKGROUND, ALWAYS };

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

enum class StereoLabels { NONE, KNOWN, ALL };

/**
 * Coloring scheme corresponding to what is available in RDKit rendering
 */
enum class ColorScheme {
    DEFAULT,
    AVALON,
    CDK,
    DARK_MODE,
    /// Black structure (to be placed against a white background)
    BLACK_WHITE,
    /// very, very light gray structure (matches DARK_MODE's carbon color)
    WHITE_BLACK
};

const std::unordered_set<ColorScheme> DARK_MODE_COLOR_SCHEMES = {
    ColorScheme::DARK_MODE, ColorScheme::WHITE_BLACK};

/**
 * The default font size for atom labels, measured in Scene pixels
 */
constexpr qreal DEFAULT_FONT_SIZE = 18;

} // namespace sketcher
} // namespace schrodinger
