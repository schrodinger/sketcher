#pragma once

#include <QtGlobal>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used to store all settings relevant to BondItems
 */
class BondItemSettings
{
  public:
    BondItemSettings() = default;
    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single BondItemSettings object with all BondItems.
    explicit BondItemSettings(const BondItemSettings&) = default;
    /// The width of the pen to use for drawing bond lines
    qreal m_bond_width = BOND_DEFAULT_PEN_WIDTH;
    /// How far apart to draw the lines of a double (or triple) bond
    qreal m_double_bond_spacing = 5.0;
    /// The width of the fat end of the wedge for chiral bonds
    qreal m_wedge_width = 5.0;
    /// The distance between hash marks in "down" (i.e. dashed) wedges
    qreal m_hash_spacing = 5.0;
};
} // namespace sketcher
} // namespace schrodinger
