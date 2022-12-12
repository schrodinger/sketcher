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
const int VIEW_SCALE = 35;

const QString FONT_NAME = "Arial";
// The default font size for atom labels (e.g. the element abbreviation),
// measured in points
const qreal DEFAULT_FONT_SIZE = 17.0;
// Ratios for the size of specified font to the size of the atom label font
const qreal SUBSCRIPT_FONT_RATIO = 0.8;
const qreal CHARGE_FONT_RATIO = 0.8;
const qreal MAPPING_FONT_RATIO = 0.6;

// The half-width of the gray rectangle drawn for predictive highlighting
const qreal PREDICTIVE_HIGHLIGHTING_HALF_WIDTH = 10;

} // namespace sketcher
} // namespace schrodinger
