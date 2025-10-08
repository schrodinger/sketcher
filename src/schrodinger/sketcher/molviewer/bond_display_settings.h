#pragma once

#include <QtGlobal>
#include <QColor>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used to store all settings relevant to BondItems
 */
class SKETCHER_API BondDisplaySettings
    : public AbstractAtomOrBondDisplaySettings
{
  public:
    BondDisplaySettings();
    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single BondDisplaySettings object with all BondItems.
    explicit BondDisplaySettings(const BondDisplaySettings&) = default;
    virtual ~BondDisplaySettings() = default;
    /// The width of the pen to use for drawing bond lines
    qreal m_bond_width = BOND_DEFAULT_PEN_WIDTH;
    /// How far apart to draw the lines of a double (or triple) bond
    qreal m_double_bond_spacing = DEFAULT_DOUBLE_BOND_SPACING;
    /// The width of the fat end of the wedge for chiral bonds
    qreal m_wedge_width = DEFAULT_BOND_WEDGE_WIDTH;
    /// The distance between hash marks in "down" (i.e. dashed) wedges
    qreal m_hash_spacing = DEFAULT_BOND_HASH_SPACING;
    /// The bond color
    QColor m_color;
    /// Whether to display stereochemistry labels when present
    bool m_stereo_labels_shown = true;

    /**
     * Scale bond width by the specified value.  Note that double bond spacing
     * and hash spacing will also be affected by this scaling.
     *
     * @param scale The scaling to apply.  A value of one will use the default
     * settings.  Values of less than one will result in narrower bonds, and
     * values greater than one will result in wider bonds.
     */
    void setScale(qreal scale);

    /**
     * Set the color scheme for the bond display settings. This affects only the
     * color of annotations (which is lighter on dark modes)
     * @param scheme The color scheme to use
     */
    void setColorScheme(const ColorScheme& scheme,
                        const QColor& carbon_color = QColor()) override;
};
} // namespace sketcher
} // namespace schrodinger
