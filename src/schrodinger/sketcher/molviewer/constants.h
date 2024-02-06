/**
 * Constants use in the molviewer Qt graphics view classes.
 */

#pragma once

#include <QColor>
#include <QString>

#include <rdkit/GraphMol/Depictor/DepictUtils.h>

#include "schrodinger/rdkit_extensions/sgroup.h"

namespace schrodinger
{
namespace sketcher
{

// the name of the rdkit property that holds user-defined colors for atoms and
// bonds
const std::string USER_COLOR = "USER_COLOR";

// The default bond length as defined by RDKit
const float BOND_LENGTH = RDDepict::BOND_LEN;

// The scaling of the graphics scene coordinate system relative to the RDKit
// coordinate system.  Without this scaling, even a 1 pt font in the graphics
// scene would be far too large for the molecule, and Qt won't render fonts
// smaller than 1 pt.
const int VIEW_SCALE = std::floor(50 / BOND_LENGTH);

const QString FONT_NAME = "Arial";
// The default font size for atom labels (e.g. the element abbreviation),
// measured in points
const qreal DEFAULT_FONT_SIZE = 13.0;
// Ratios for the size of specified font to the size of the atom label font
const qreal SUBSCRIPT_FONT_RATIO = 0.6;
const qreal CHARGE_FONT_RATIO = 0.6;
const qreal MAPPING_FONT_RATIO = 0.5;
const qreal CHIRALITY_FONT_RATIO = 0.6;
const qreal SGROUP_FONT_RATIO = 1.0;
const qreal CURSOR_HINT_FONT_RATIO = 1.0;

// The font size used to display the number of bonds that the Draw Chain scene
// tool will draw
const qreal DRAW_CHAIN_FONT_SIZE = 11.0;

// Ratios between the dot width for radicals and the font size
const qreal RADICAL_DOT_RATIO = 0.25;

// margin around the text in an atom text rectangle in pixels
const qreal TEXT_RECT_MARGIN = 1.3;

// the hight ratio on the H label at which the top of the h count label is drawn
const qreal H_COUNT_LABEL_VERTICAL_OFFSET_RATIO = 0.7;

// Ratio between the distance of radical dots from the center and the label
// height. Used for radicals when the atom label is visible
const qreal RADICAL_DISTANCE_FROM_LABEL_RATIO = 0.75;
// Ratio between the distance of radical dots from the center and the radical
// dot size. Used for radicals when the atom label is not visible
const qreal RADICAL_DISTANCE_FROM_HIDDEN_LABEL_RATIO = 3.;

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

// the color of structure hints (e.g. the blue lines that show where a bond will
// be drawn when using the bond tool or the fragment structure when mousing over
// a bond when using the ring/fragment tool)
const QColor STRUCTURE_HINT_COLOR = QColor("#9cbcd1");

// the color of the cursor hint (e.g. the element name that attaches to the
// mouse cursor when using the atom tool, or the bond icon that attaches to the
// cursor when using the bond tool)
const QColor CURSOR_HINT_COLOR = QColor("#5191bb");

// This should be set to the shade of gray used in the tool button icons.
// cursor_hint_from_svg in scene_utils will replace this color with
// CURSOR_HINT_COLOR.
const QRgb TOOL_BUTTON_DARK_GRAY = qRgb(102, 102, 102);

// The width of bonds for fragment hints (e.g. the blue rings that appear when
// using the ring tool)
const qreal FRAGMENT_HINT_BOND_WIDTH = 1.0;

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

// the ratio of the scroll step for scrolling with the key arrows to a bond
// length
const float KEY_SCROLL_BOND_LENGTH_RATIO = 0.5;

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

// The maximum and minimum charg that can be set on an atom
const int MAX_ATOM_CHARGE = 8;
const int MIN_ATOM_CHARGE = -8;

const QColor ROTATION_ITEM_COLOR = QColor("#ff9b00");

const float ROTATION_ITEM_HANDLE_RADIUS = 12.f;
const float ROTATION_ITEM_PIVOT_RADIUS = 8.f;
const float ROTATION_ITEM_PEN_WIDTH = 3.f;

// The font size and color used to display the angle value for the rotation item
const qreal ROTATION_ITEM_FONT_SIZE = 11.0;
const QColor ROTATION_ITEM_TEXT_COLOR = CHIRALITY_LABEL_COLOR;

// The width of the pen used for painting non-molecular objects (reaction arrows
// and pluses) measured in Scene units
const qreal NON_MOLECULAR_PEN_WIDTH = 3.0;

// The dimensions of the reaction arrow, measured in Scene units
const qreal ARROW_LENGTH = 40.0;
const qreal ARROW_WIDTH = 10.0;

// The dimensions of the reaction plus sign, measured in Scene units
const qreal PLUS_LENGTH = 20.0;

// The *_SPACING values indicate, when importing a molecule or a reaction, how
// much space to leave between...  (measured in RDKit units)

// ...two reactants or two products of a reaction
const double PLUS_SPACING = BOND_LENGTH * 1.5;

// ...the last reactant and the first product of a reaction
const double ARROW_SPACING = PLUS_SPACING * 2;

// ...the current model contents and the newly imported molecule (of first
// molecule of the reaction)
const double IMPORT_SPACING = BOND_LENGTH * 2;

// The dimensions of the atom mapping tool arrow, measured in Scene units
const qreal MAPPING_ARROW_HEAD_LENGTH = 10;
const qreal MAPPING_ARROW_HEAD_HALF_WIDTH = 3;

// The width of the path drawn around non-molecular objects for predictive and
// selection highlighting
const qreal NON_MOLECULAR_HIGHLIGHT_PADDING = 6.0;

// The width of the path drawn around substance groups for predictive and
// selection highlighting
const qreal S_GROUP_HIGHLIGHT_PADDING = 6.0;

// When inserting a fragment (e.g. the Draw Ring tool), we look for fragment
// atoms that have identical coordinates to core atoms.  This is the max
// distance (in MolModel coordinates) for two coordinates to be considered
// "identical."
const double MAX_DIST_FOR_FRAG_OVERLAP = 0.01;

// the dimensions of brackets around substance groups. BRAKETS_LONG_SIDE is
// defined in rdkit_extensions/sgroup.h
const qreal BRACKETS_SHORT_SIDE = rdkit_extensions::BRACKETS_LONG_SIDE * 0.2;
const qreal LABEL_DISTANCE_FROM_BRACKETS = BRACKETS_SHORT_SIDE * 0.5;

// At the end of a drag-and-drop with the Move/Rotate tool, overlapping atoms
// from different molecules will be merged.  This is the maximum distance (and
// squared maximum distance) for two atoms to be considered overlapping.  Both
// values are in MolModel coordinates.
const double MAX_DIST_FOR_DRAG_MERGE = 0.2 * BOND_LENGTH;
const double MAX_DIST_SQ_FOR_DRAG_MERGE =
    MAX_DIST_FOR_DRAG_MERGE * MAX_DIST_FOR_DRAG_MERGE;

// The radius of the circle drawn around overlapping atoms during a
// drag-and-drop with the Merge/Rotate tool.  Given in Scene coordinates.
const qreal DRAG_MERGE_HINT_RADIUS = VIEW_SCALE * MAX_DIST_FOR_DRAG_MERGE;

// The thickness of the pen used to draw the circles drawn around overlapping
// atoms during a drag-and-drop with the Merge/Rotate tool.  Given in Scene
// coordinates.
const qreal DRAG_MERGE_HINT_WIDTH = 4.0;

// The image to use for the standard arrow mouse cursor
const QString ARROW_CURSOR_PATH = ":/icons/cursor-arrow.svg";

// The "hotspot" of the above cursor image, i.e., (x, y) coordinates for where
// in the image the click should happen
const int CURSOR_HOTSPOT_X = 2;
const int CURSOR_HOTSPOT_Y = 5;

// When attaching a hint image to the cursor, the (x, y) coordinates of the
// upper-left corner of the hint
const int CURSOR_HINT_X = 16;
const int CURSOR_HINT_Y = 18;

// The maximum allowable width and height for a cursor hint image generated from
// a toolbar icon.  Note that the cursor hint image won't necessarily be a
// square of this exact size due to both the aspect ratio and because the image
// is cropped after resizing.
const int CURSOR_HINT_IMAGE_SIZE = 28;

// The Z ordering for graphics items.  Items listed later in this enum (i.e. a
// higher value) will be drawn on top of items listed earlier in this enum (i.e.
// a lower value)
enum class ZOrder {
    SELECTION_HIGHLIGHTING = 1,
    PREDICTIVE_HIGHLIGHTING,
    BOND,
    ATOM,
    SGROUP,
    RXN_ARROW_AND_PLUS,
    DRAG_SELECT_OUTLINE,
    HINT,
    ROTATION_HANDLE,
};

// The atomic number of RDKit dummy atoms
const int DUMMY_ATOMIC_NUMBER = 0;

// Bit flags for specifying subsets of Scene items based on the type of model
// object they represent
typedef uint8_t InteractiveItemFlagType;
namespace InteractiveItemFlag
{
enum : InteractiveItemFlagType { // clang-format off
    ATOM_NOT_R_NOT_AP     = 1 << 0,
    R_GROUP               = 1 << 1,
    S_GROUP               = 1 << 2,
    ATTACHMENT_POINT      = 1 << 3,
    BOND_NOT_AP           = 1 << 4,
    ATTACHMENT_POINT_BOND = 1 << 5,
    NON_MOLECULAR         = 1 << 6, // clang-format on
    ATOM_NOT_AP = ATOM_NOT_R_NOT_AP | R_GROUP,
    ATOM = ATOM_NOT_AP | ATTACHMENT_POINT,
    BOND = BOND_NOT_AP | ATTACHMENT_POINT_BOND,
    MOLECULAR = ATOM | BOND | S_GROUP,
    ALL = MOLECULAR | NON_MOLECULAR,
    ATTACHMENT_POINT_OR_AP_BOND = ATTACHMENT_POINT | ATTACHMENT_POINT_BOND,
    MOLECULAR_NOT_AP = ATOM_NOT_AP | BOND_NOT_AP,
};
}

} // namespace sketcher
} // namespace schrodinger
