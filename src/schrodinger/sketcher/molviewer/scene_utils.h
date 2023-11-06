#pragma once

#include <vector>
#include <unordered_map>

#include <QList>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/constants.h"

class QColor;
class QFont;
class QGraphicsItem;
class QPixmap;
class QString;

namespace RDKit
{
class Atom;
class Bond;
class Conformer;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class AtomItem;
class AtomItemSettings;
class BondItem;
class BondItemSettings;
class Bonds;
class Fonts;

/**
 * Create all graphics items needed to represent the given molecule
 * @param mol The molecule to draw
 * @param fonts The fonts to use for rendering the molecule
 * @param atom_item_settings The settings for displaying atoms
 * @param bond_item_settings The settings for displaying bonds
 * @param draw_attachment_points Whether to create graphics items for
 * attachment point atoms.
 * @return A tuple of
 *   - A list of all newly created graphics items.
 *   - A map of atom -> the graphics item used to represent that atom
 *   - A map of bond -> the graphics item used to represent that bond
 */
SKETCHER_API std::tuple<std::vector<QGraphicsItem*>,
                        std::unordered_map<const RDKit::Atom*, AtomItem*>,
                        std::unordered_map<const RDKit::Bond*, BondItem*>>
create_graphics_items_for_mol(const RDKit::ROMol* mol, const Fonts& fonts,
                              AtomItemSettings& atom_item_settings,
                              BondItemSettings& bond_item_settings,
                              const bool draw_attachment_points = true);

/**
 * Update all graphics items to represent an updated conformer
 * @param atom_items All atom items to update
 * @param bond_items All bond items to update
 * @param mol The molecule containing the updated conformer
 */
SKETCHER_API void
update_conf_for_mol_graphics_items(const QList<QGraphicsItem*>& atom_items,
                                   const QList<QGraphicsItem*>& bond_items,
                                   const RDKit::ROMol& mol);

/**
 * Draw the given text to an appropriately sized pixmap
 * @param text The text to draw
 * @param font The font to use
 * @param color The desired color of the text
 * @return The newly generated pixmap
 */
SKETCHER_API QPixmap render_text_to_pixmap(const QString& text,
                                           const QFont& font,
                                           const QColor& color);

/**
 * Create a pixmap for use as a cursor hint from the given tool button icon.
 * The returned pixmap will be sized correctly and (optionally) recolored to the
 * appropriate shade of blue.
 * @param path The path to the SVG tool button image.  This can be either a
 * filesystem path or a Qt resource path.
 * @param recolor Whether to convert the all TOOL_BUTTON_DARK_GRAY pixels to
 * CURSOR_HINT_COLOR
 * @return The newly generated pixmap
 */
SKETCHER_API QPixmap cursor_hint_from_svg(const QString& path,
                                          const bool recolor = true);

/**
 * @return A pixmap containing the arrow cursor
 */
SKETCHER_API QPixmap get_arrow_cursor_pixmap();

} // namespace sketcher
} // namespace schrodinger
