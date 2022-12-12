/**
 * Constants use in the molviewer Qt graphics view classes.
 *
 * Copyright Schrodinger, LLC. All rights reserved.
 */

#pragma once

#include <QGraphicsItem>
#include <QString>

namespace schrodinger
{
namespace sketcher
{

// The Angstrom to point (i.e. font size) conversion factor.  Coordgen (via
// RDKit) attempts to set bond lengths to a reasonable number of Angstroms, even
// in 2D.  As a result, without this scaling, the default font size would be
// substantially larger than an entire molecule and even a 1 pt font would be
// far too large.
const int VIEW_SCALE = 50;

const QString FONT_NAME = "Arial";
// The default font size for atom labels (e.g. the element abbreviation),
// measured in points
const qreal DEFAULT_FONT_SIZE = 13.0;
// Ratios for the size of specified font to the size of the atom label font
const qreal SUBSCRIPT_FONT_RATIO = 0.8;
const qreal CHARGE_FONT_RATIO = 0.8;
const qreal MAPPING_FONT_RATIO = 0.6;

// The half-width of the gray rectangle drawn for predictive highlighting
const qreal PREDICTIVE_HIGHLIGHTING_HALF_WIDTH = 10;

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

} // namespace sketcher
} // namespace schrodinger
