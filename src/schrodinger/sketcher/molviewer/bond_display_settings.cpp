#include "schrodinger/sketcher/molviewer/bond_display_settings.h"

#include <stdexcept>
#include <stdio.h>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

qreal static scale_value(qreal scale, qreal default_value, qreal scaling_factor)
{
    return (1 + scaling_factor * (scale - 1)) * default_value;
}

void BondDisplaySettings::setScale(qreal scale)
{
    if (scale <= 0) {
        throw std::invalid_argument("Bond width scale must be positive");
    }
    m_bond_width = scale * BOND_DEFAULT_PEN_WIDTH;
    m_double_bond_spacing = scale_value(scale, DEFAULT_DOUBLE_BOND_SPACING,
                                        DOUBLE_BOND_SPACING_SCALING_FACTOR);
    m_hash_spacing = scale_value(scale, DEFAULT_BOND_HASH_SPACING,
                                 BOND_HASH_SPACING_SCALING_FACTOR);
}

} // namespace sketcher
} // namespace schrodinger
