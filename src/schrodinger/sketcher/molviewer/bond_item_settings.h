/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <QtGlobal>

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
    qreal m_bond_width = 2.2;
    /// How far apart to draw the lines of a double (or triple) bond
    qreal m_double_bond_spacing = 5;
};
} // namespace sketcher
} // namespace schrodinger
