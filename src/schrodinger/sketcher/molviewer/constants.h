/**
 * Constants use in the molviewer Qt graphics view classes.
 *
 * Copyright Schrodinger, LLC. All rights reserved.
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
const qreal SUBSCRIPT_FONT_RATIO = 0.8;
const qreal CHARGE_FONT_RATIO = 0.8;
const qreal MAPPING_FONT_RATIO = 0.6;

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

// The Z ordering for graphics items.  Items listed later in this enum (i.e. a
// higher value) will be drawn on top of items listed earlier in this enum (i.e.
// a lower value)
enum class ZOrder {
    SELECTION_HIGHLIGHTING = 1,
    PREDICTIVE_HIGHLIGHTING,
    BOND,
    ATOM,
    DRAG_SELECT_OUTLINE,
};

} // namespace sketcher
} // namespace schrodinger
