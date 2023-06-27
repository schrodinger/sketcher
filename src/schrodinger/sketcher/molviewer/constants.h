/**
 * Constants use in the molviewer Qt graphics view classes.
 *
 */

#pragma once

#include <QColor>
#include <QString>

namespace schrodinger
{
namespace sketcher
{

// The scaling of the graphics scene coordinate system relative to the RDKit
// coordinate system.  Without this scaling, even a 1 pt font in the graphics
// scene would be far too large for the molecule, and Qt won't render fonts
// smaller than 1 pt.
const int VIEW_SCALE = 50;

const QString FONT_NAME = "Arial";
// The default font size for atom labels (e.g. the element abbreviation),
// measured in points
const qreal DEFAULT_FONT_SIZE = 13.0;
// Ratios for the size of specified font to the size of the atom label font
const qreal SUBSCRIPT_FONT_RATIO = 0.6;
const qreal CHARGE_FONT_RATIO = 0.6;
const qreal MAPPING_FONT_RATIO = 0.6;
const qreal CHIRALITY_FONT_RATIO = 0.6;

// Ratios between the dot width for radicals and the font size
const qreal RADICAL_DOT_RATIO = 0.25;

// Ratios the spacer between two consecutive labels (such as NH) and the width
// of the main element label
const qreal LABEL_SPACER_RATIO = 0.1;

// Ratio between the distance of radical dots from the center and the label
// height. Used for radicals when the atom label is visible
const qreal RADICAL_DISTANCE_FROM_LABEL_RATIO = 0.75;
// Ratio between the distance of radical dots from the center and the radical
// dot size. Used for radicals when the atom label is not visible
const qreal RADICAL_DISTANCE_FROM_HIDDEN_LABEL_RATIO = 3.;

// QFontmetrics returns a bigger height for characters that the actual size.
// This factor is used to correct that so the label fits the actual dimensions
// of the text
const qreal LABEL_RECT_HEIGHT_ADJUSTMENT_FACTOR = 0.2;

// Ratios between paired dots and their size for radicals
const qreal PAIRED_ELECTRONS_DISTANCE_RATIO = 0.75;

// The radius of the circle drawn around atoms for predictive highlighting
const qreal ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS = 10;

// The half-width of the rectangle drawn around bonds for predictive
// highlighting
const qreal BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH = 10;

// The radius of the circle drawn around atoms for selection highlighting
const qreal ATOM_SELECTION_HIGHLIGHTING_RADIUS = 11;

// The half-width of the rectangle drawn around bonds for selection highlighting
const qreal BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH = 8;

// The dimensions of the arrow for dative bonds (a.k.a. coordinate bonds)
const qreal DATIVE_ARROW_HALF_WIDTH = 4.0;
const qreal DATIVE_ARROW_LENGTH = 7.0;

// The width of the margin around atom labels, i.e. how many pixels should we
// leave between the end of the bond and the start of the atom label
const qreal ATOM_LABEL_MARGIN = 4.0;

// By how many pixels should we shorten each end of the inner line of a double
// or aromatic bond (e.g. the line that's drawn on the inside of the ring)
const qreal DOUBLE_BOND_INNER_LINE_SHORTENING = 6.0;

// When picking the "best ring" for a double bond or an aromatic bond (i.e. the
// ring that the second line (or the dashed line) of the bond should be drawn
// inside of), rings with DOUBLE_BOND_BEST_RING_SIZE_CUTOFF atoms or fewer will
// always be preferred over rings with more than
// DOUBLE_BOND_BEST_RING_SIZE_CUTOFF atoms.
const unsigned int DOUBLE_BOND_BEST_RING_SIZE_CUTOFF = 8;

// The background color for drawing selection highlighting
const QColor SELECTION_BACKGROUND_COLOR = QColor("#c7d5b8");

// The outline color for drawing selection highlighting
const QColor SELECTION_OUTLINE_COLOR = QColor("#779c59");

// The color for drawing predictive highlighting.  The same color is used for
// both the outline and the background.
const QColor PREDICTIVE_HIGHLIGHTING_COLOR = QColor("#e2eadb");

// the color of the dotted line around valence errors
const QColor VALENCE_ERROR_BORDER_COLOR = QColor("#fb7100");

// the color of the stereo annotations (e.g. "R" and "S")
const QColor CHIRALITY_LABEL_COLOR = QColor("#333333");

// the color of hint bonds (e.g. the blue lines that show where a bond will be
// drawn)
const QColor HINT_COLOR = Qt::blue;

// the distance between the stereo annotation and the atom in bond length units
const qreal CHIRALITY_LABEL_DISTANCE_RATIO = 0.25;

// the distance between the stereo annotation and the bond in bond length units
const qreal BOND_STEREO_LABEL_DISTANCE_RATIO = 0.20;

// the size of the dotted line around valence errors
const qreal VALENCE_ERROR_BORDER_WIDTH = 2.f;

// the color of the area of valence errors
const QColor VALENCE_ERROR_AREA_COLOR = QColor("#ffecc5");

// the border size around the main element to define the valence error area
const qreal VALENCE_ERROR_AREA_BORDER = 3.f;

// the amount the scene should be zoomed  out to allow for some white space
// around a molecule when they fill the screen
const float FIT_TO_SCREEN_MARGIN_FACTOR = 0.15f;

/// The width of the pen to use for drawing bond lines
const qreal BOND_DEFAULT_PEN_WIDTH = 2.2;

// The length of the bond to an attachment point in bond length units
const qreal ATTACHMENT_POINT_BOND_DISTANCE_RATIO = 0.5;

// The number of "waves" in the squiggle used to represent attachment points
// (where each wave is goes from zero, to the maximum, to the minimum, and back
// to zero)
const int ATTACHMENT_POINT_SQUIGGLE_NUMBER_OF_WAVES = 3;

// The width of each "wave" in the squiggle used to represent attachment points.
// This value is measured in Scene units.
const qreal ATTACHMENT_POINT_SQUIGGLE_WIDTH_PER_WAVE = 8.0;

// The height of the attachment point squiggle, measured in Scene units
const qreal ATTACHMENT_POINT_SQUIGGLE_HEIGHT = 3.0;

// The maximum number of attachment points that can be addded to any given atom
// when using the attachment point tool.  Note, however, that it's possible to
// *import* structures that violate this limit.
const unsigned int MAX_ATTACHMENT_POINTS_PER_ATOM = 2;

// The Z ordering for graphics items.  Items listed later in this enum (i.e. a
// higher value) will be drawn on top of items listed earlier in this enum (i.e.
// a lower value)
enum class ZOrder {
    SELECTION_HIGHLIGHTING = 1,
    PREDICTIVE_HIGHLIGHTING,
    BOND,
    ATOM,
    DRAG_SELECT_OUTLINE,
    HINT,
};

// The atomic number of RDKit dummy atoms
const int DUMMY_ATOMIC_NUMBER = 0;

/**
 * This is the reference bond length we use when rescaling input molecules.
 * The current value of 1.17 comes from playing a bit with coordgen.
 */
const double DEFAULT_MOLVIEWER_BOND_LENGTH = 1.17;

} // namespace sketcher
} // namespace schrodinger
