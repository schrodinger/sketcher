/**
 * Constants use in the molviewer Qt graphics view classes.
 */

#pragma once

#include <QColor>
#include <QString>

#include <rdkit/GraphMol/Depictor/DepictUtils.h>

#include "schrodinger/sketcher/rdkit/sgroup.h"

namespace schrodinger
{
namespace sketcher
{

/**
 *the color of stereo annotations on chiral atoms (e.g. (r) or (s) )
 */
const auto STEREO_ANNOTATION_COLOR = QColor("#333333");

/**
 * Maximum and minimum ammount of charge that
 * can be assigned to an atom
 */
constexpr int ATOM_CHARGE_LIMIT = 8;

/**
 * Max value for Isotope
 * that can be assigned to an atom
 */
constexpr int MAX_ISOTOPE_VALUE = 999;

/**
 * Maximum and minimum numbers of unpaired electrons that can be assigned to an
 * atom
 */
constexpr int MIN_UNPAIRED_E = 0;
constexpr int MAX_UNPAIRED_E = 4;

// the name of the rdkit property that holds user-defined colors for atoms and
// bonds
const std::string USER_COLOR = "USER_COLOR";

// The default bond length as defined by RDKit
const float BOND_LENGTH = RDDepict::BOND_LEN;

// the width of empty space around the scene items in a saved picture (in
// pixels)
const int SAVED_PICTURE_PADDING = 10;

// The scaling of the graphics scene coordinate system relative to the RDKit
// coordinate system.  Without this scaling, even a 1 pt font in the graphics
// scene would be far too large for the molecule, and Qt won't render fonts
// smaller than 1 pt.
const int VIEW_SCALE = std::floor(50 / BOND_LENGTH);

// Mac uses a smaller default system cursor than other platforms, so we reduce
// the size of our cursor on Mac
#ifdef __APPLE__
const qreal CURSOR_SCALE = 0.8;
#else
const qreal CURSOR_SCALE = 1.0;
#endif

const QString FONT_NAME = "Arimo";
// Ratios for the size of specified font to the size of the atom label font.
// Note that the atom label font size is declared in the public image_constants
// header (mmshare/include/schrodinger/sketcher/image_constants.h) because it
// needs to be accessed from image_generation.h (which is also public)
const qreal SUBSCRIPT_FONT_RATIO = 0.6;
const qreal CHARGE_FONT_RATIO = 0.6;
const qreal MAPPING_FONT_RATIO = 0.5;
const qreal QUERY_LABEL_FONT_RATIO = 0.5;
const qreal CHIRALITY_FONT_RATIO = 0.6;
const qreal SGROUP_FONT_RATIO = 1.0;
const qreal CURSOR_HINT_FONT_RATIO = 1.0 * CURSOR_SCALE;

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

// the distance between an atom and its query label in bond length units
const qreal QUERY_LABEL_DISTANCE_RATIO = 0.35;

// the maximum length of characters displayed for a SMARTS label in scene
// coordinates. If the label is longer than this, it will be truncated and an
// ellipsis will be added
const qreal MAX_SMARTS_LABEL_LENGTH = 1.2 * BOND_LENGTH * VIEW_SCALE;

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

const QColor SELECT_TOOL_LINE_COLOR_DARK_BG = QColor("#eeeeee");
const QColor SELECT_TOOL_LINE_COLOR = QColor("#000000");

const QColor DARK_BACKGROUND_COLOR = QColor("#1e1e1e");
const QColor LIGHT_BACKGROUND_COLOR = QColor("#ffffff");

// The background color for drawing selection highlighting
const QColor SELECTION_BACKGROUND_COLOR = QColor("#c7d5b8");
const QColor SELECTION_BACKGROUND_COLOR_DARK_BG = QColor("#4d5d3f");

// The outline color for drawing selection highlighting
const QColor SELECTION_OUTLINE_COLOR = QColor("#779c59");
const QColor SELECTION_OUTLINE_COLOR_DARK_BG = QColor("#82a664");

// The color for drawing predictive highlighting.  The same color is used for
// both the outline and the background.
const QColor PREDICTIVE_HIGHLIGHTING_COLOR = QColor("#e2eadb");
const QColor PREDICTIVE_HIGHLIGHTING_COLOR_DARK_BG = QColor("#555a54");

// the color for drawing predictive highlighting for the move-rotate tool
const QColor MOVE_ROTATE_PREDICTIVE_HIGHLIGHTING_COLOR = QColor("#eaeaea");
const QColor MOVE_ROTATE_PREDICTIVE_HIGHLIGHTING_COLOR_DARK_BG =
    QColor("#4c4c4c");

// the color of the dotted line around valence errors
const QColor VALENCE_ERROR_BORDER_COLOR = QColor("#fb7100");

// the color of the stereo annotations (e.g. "R" and "S")
const QColor ANNOTATION_COLOR = QColor("#333333");
const QColor ANNOTATION_COLOR_DARK_BG = QColor("#bababa");

// the color of structure hints (e.g. the blue lines that show where a bond will
// be drawn when using the bond tool or the fragment structure when mousing over
// a bond when using the ring/fragment tool)
const QColor STRUCTURE_HINT_COLOR = QColor("#9cbcd1");

// the color of the cursor hint (e.g. the element name that attaches to the
// mouse cursor when using the atom tool, or the bond icon that attaches to the
// cursor when using the bond tool)
const QColor CURSOR_HINT_COLOR = QColor("#5191bb");

// the color of the line that connects an allowed/disallowed query list or
// SMARTS query to the atom it refers to
const QColor QUERY_LABEL_LINE_COLOR = QColor(102, 102, 102);

// the color of query bonds
const QColor QUERY_BOND_COLOR = QUERY_LABEL_LINE_COLOR;

// This should be set to the shade of gray used in the tool button icons.
// cursor_hint_from_svg in scene_utils will replace this color with
// CURSOR_HINT_COLOR.
const QRgb TOOL_BUTTON_DARK_GRAY = qRgb(102, 102, 102);

// The width of bonds for fragment hints (e.g. the blue rings that appear when
// using the ring tool and hovering over an atom or bond)
const qreal FRAGMENT_HINT_BOND_WIDTH = 1.0;

// The width of bonds and the bond line spacing for the fragment cursor pixmap
// (i.e. the small fragment image that follows the mouse cursor when hovering
// over empty space).  Note that these values should be much larger than the
// typical bond width and spacing because the image will be scaled down to only
// CURSOR_HINT_IMAGE_SIZE x CURSOR_HINT_IMAGE_SIZE pixels.
const qreal FRAGMENT_CURSOR_HINT_BOND_WIDTH = 7.0;
const qreal FRAGMENT_CURSOR_HINT_BOND_SPACING = 12.0;

// the distance between the (circle circumscribed to the) stereo annotation and
// the atom in bond length units.
const qreal CHIRALITY_LABEL_DISTANCE_RATIO = 0.10;

// the distance between the stereo annotation and the bond in bond length units
const qreal BOND_STEREO_LABEL_DISTANCE_RATIO = 0.20;

// the size of the dotted line around valence errors
const qreal VALENCE_ERROR_BORDER_WIDTH = 2.f;

// the color of the area of valence errors
const QColor VALENCE_ERROR_AREA_COLOR = QColor("#ffecc5");
const QColor VALENCE_ERROR_AREA_COLOR_DARK_BG = QColor("#71471e");

// the border size around the main element to define the valence error area
const qreal VALENCE_ERROR_AREA_BORDER = 3.f;

// the amount the scene should be zoomed  out to allow for some white space
// around a molecule when they fill the screen
const float FIT_TO_SCREEN_MARGIN_FACTOR = 0.15f;

// the ratio of the scroll step for scrolling with the key arrows to a bond
// length
const float KEY_SCROLL_BOND_LENGTH_RATIO = 0.5;

/// The width of the pen to use for drawing bond lines
const qreal BOND_DEFAULT_PEN_WIDTH = 2.4;

/// The minimum length of a bond, measured in Scene units. Avoiding collisions
/// with the atom labels will not force bonds to be shorter than this.
const qreal MINIMUM_BOND_LENGTH = BOND_LENGTH * VIEW_SCALE * 0.10;

/// The default spacing between the two lines of a double bond, or between pairs
/// of adjacent lines in triple bonds.
const qreal DEFAULT_DOUBLE_BOND_SPACING = 5.5;

/// The default width of the fat end of the wedge for wedge (up or down) bonds
const qreal DEFAULT_BOND_WEDGE_WIDTH = 6.0;

// The default spacing between hash marks in a down wedge bond
const qreal DEFAULT_BOND_HASH_SPACING = 5.0;

/// How double bond spacing is scaled relative to bond width. If this value is
/// one, then double bond spacing will be scaled the same as bond width. E.g.,
/// if bond width is doubled, then double bond spacing will be doubled. If this
/// value is less than one, then changes in bond width will be de-emphasized for
/// double bond spacing. E.g., doubling the bond width will still increase the
/// double bond spacing, but by _less_ than a doubling. If this value is
/// greater than one, then changes will be exaggerated (so _more_ than a
/// doubling in our example). The greater the difference between this value and
/// one, the more pronounced this effect will be.
///
/// Note that this value must be positive.
const qreal DOUBLE_BOND_SPACING_SCALING_FACTOR = 0.9;

/// How hash spacing (for down wedge bonds) is scaled relative to the bond
/// width. See DOUBLE_BOND_SPACING_SCALING_FACTOR for an explanation of how this
/// scaling works.
const qreal BOND_HASH_SPACING_SCALING_FACTOR = 0.3;

/// The opacity of the portion of bonds that are behind bond labels
const qreal OPACITY_OF_BOND_BEHIND_LABEL = 0.15;

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

const QColor ROTATION_ITEM_COLOR = QColor("#ff9b00");

const float ROTATION_ITEM_HANDLE_RADIUS = 12.f;
const float ROTATION_ITEM_PIVOT_RADIUS = 8.f;
const float ROTATION_ITEM_PEN_WIDTH = 3.f;

// The font size and color used to display the angle value for the rotation item
const qreal ROTATION_ITEM_FONT_SIZE = 11.0;
const QColor ROTATION_ITEM_TEXT_COLOR = ANNOTATION_COLOR;

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
// defined in sketcher/rdkit/sgroup.h
const qreal BRACKETS_SHORT_SIDE = BRACKETS_LONG_SIDE * 0.2;
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
// in the image the click should happen. (Note that we don't bother to multiply
// these values by CURSOR_SCALE because we can't put 0.8 into an integer.
// Instead, we translate the image in get_arrow_cursor_pixmap to keep the
// hotspot in the right place.)
const int CURSOR_HOTSPOT_X = 1;
const int CURSOR_HOTSPOT_Y = 1;

// When attaching a hint image to the cursor, the (x, y) coordinates of the
// upper-left corner of the hint
const int CURSOR_HINT_X = 15 * CURSOR_SCALE;
const int CURSOR_HINT_Y = 13 * CURSOR_SCALE;

// The maximum allowable width and height for a cursor hint image generated from
// a toolbar icon.  Note that the cursor hint image won't necessarily be a
// square of this exact size due to both the aspect ratio and because the image
// is cropped after resizing.
const int CURSOR_HINT_IMAGE_SIZE = 28 * CURSOR_SCALE;

// the colors for the cursor arrow
const QColor LIGHT_CURSOR_COLOR = DARK_BACKGROUND_COLOR;
const QColor DARK_CURSOR_COLOR = LIGHT_BACKGROUND_COLOR;

// The maximum number of atoms that can be imported at once. Attempting to
// import a molecule or reaction with too many atoms will result in an error
// dialog informing the user of the problem.
const unsigned int MAX_NUM_ATOMS_FOR_IMPORT = 200;

// maximum allowed values for spin boxes in the Edit Aom Properties dialog
const unsigned int MAX_STEREO_GROUP_ID = 99;
const unsigned int MAX_QUERY_TOTAL_H = 99;
const unsigned int MAX_QUERY_NUM_CONNECTIONS = 99;
const unsigned int MAX_QUERY_RING_COUNT = 99;
const unsigned int MAX_QUERY_RING_BOND_COUNT = 99;
const unsigned int MAX_QUERY_SMALLEST_RING_SIZE = 999;

// the icon to use for the clear button on line edits and BlankableSpinBoxes
const QString LINE_EDIT_CLEAR_ICON_PATH = ":icons/clear_line_edit.png";

// The Z ordering for graphics items.  Items listed later in this enum (i.e. a
// higher value) will be drawn on top of items listed earlier in this enum (i.e.
// a lower value)
enum class ZOrder {
    BOND_HIGHLIGHTING = 1,
    ATOM_HIGHLIGHTING,
    SELECTION_HIGHLIGHTING,
    PREDICTIVE_HIGHLIGHTING,
    MONOMER_CONNECTOR,
    BOND,
    MONOMER,
    ATOM,
    SGROUP,
    RXN_ARROW_AND_PLUS,
    DRAG_SELECT_OUTLINE,
    HINT,
    ROTATION_HANDLE,
};

// Bit flags for specifying subsets of Scene items based on the type of model
// object they represent
typedef uint16_t InteractiveItemFlagType;
namespace InteractiveItemFlag
{
enum : InteractiveItemFlagType { // clang-format off
    NONE                  = 0,
    ATOM_NOT_R_NOT_AP     = 1 << 0,
    R_GROUP               = 1 << 1,
    S_GROUP               = 1 << 2,
    ATTACHMENT_POINT      = 1 << 3,
    AA_MONOMER            = 1 << 4,
    NA_PHOSPHATE          = 1 << 5,
    NA_SUGAR              = 1 << 6,
    NA_BASE               = 1 << 7,
    CHEM_MONOMER          = 1 << 8,
    BOND_NOT_AP           = 1 << 9,
    ATTACHMENT_POINT_BOND = 1 << 10,
    MONOMER_CONNECTOR     = 1 << 11,
    NON_MOLECULAR         = 1 << 12, // clang-format on
    ATOM_NOT_AP = ATOM_NOT_R_NOT_AP | R_GROUP,
    ATOM = ATOM_NOT_AP | ATTACHMENT_POINT,
    BOND = BOND_NOT_AP | ATTACHMENT_POINT_BOND,
    MOLECULAR = ATOM | BOND | S_GROUP,
    ATTACHMENT_POINT_OR_AP_BOND = ATTACHMENT_POINT | ATTACHMENT_POINT_BOND,
    MOLECULAR_NOT_AP = ATOM_NOT_AP | BOND_NOT_AP,
    NA_MONOMER = NA_PHOSPHATE | NA_SUGAR | NA_BASE,
    MONOMER = AA_MONOMER | NA_MONOMER | CHEM_MONOMER,
    MONOMERIC = MONOMER | MONOMER_CONNECTOR,
    MOLECULAR_OR_MONOMERIC = MOLECULAR | MONOMERIC,
    ATOM_OR_MONOMER = ATOM | MONOMER,
    BOND_OR_CONNECTOR = BOND | MONOMER_CONNECTOR,
    ALL = MOLECULAR_OR_MONOMERIC | NON_MOLECULAR,
};
}

} // namespace sketcher
} // namespace schrodinger
