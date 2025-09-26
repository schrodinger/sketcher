#pragma once

#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <QtGlobal>
#include <QList>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/constants.h"

class QColor;
class QFont;
class QGraphicsItem;
class QLineF;
class QPainterPath;
class QPixmap;
class QString;

namespace RDKit
{
class Atom;
class Bond;
class Conformer;
class ROMol;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class AbstractAtomOrMonomerItem;
class AbstractGraphicsItem;
class AtomItem;
class AtomDisplaySettings;
class BondItem;
class BondDisplaySettings;
class Fonts;
class SGroupItem;
class NonMolecularObject;

/**
 * Types of monomers
 */
enum class MonomerType { PEPTIDE, NA_BASE, NA_PHOSPHATE, NA_SUGAR, CHEM };

/**
 * Determine what type of monomer the given atom represents.
 *
 * @throw std::runtime_error if the atom does not represent a monomer
 */
MonomerType get_monomer_type(const RDKit::Atom* atom);

/**
 * Create all graphics items needed to represent the given molecule
 * @param mol The molecule to draw
 * @param fonts The fonts to use for rendering the molecule
 * @param atom_display_settings The settings for displaying atoms
 * @param bond_display_settings The settings for displaying bonds
 * @param draw_attachment_points Whether to create graphics items for
 * attachment point atoms.
 * @return A tuple of
 *   - A list of all newly created graphics items.
 *   - A map of atom -> the graphics item used to represent that atom
 *   - A map of bond -> the graphics item used to represent that bond
 *   - A map of S-group -> the graphics item used to represent that S-group
 */
SKETCHER_API
std::tuple<std::vector<QGraphicsItem*>,
           std::unordered_map<const RDKit::Atom*, QGraphicsItem*>,
           std::unordered_map<const RDKit::Bond*, QGraphicsItem*>,
           std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>>
create_graphics_items_for_mol(const RDKit::ROMol* mol, const Fonts& fonts,
                              const AtomDisplaySettings& atom_display_settings,
                              const BondDisplaySettings& bond_display_settings,
                              const bool draw_attachment_points = true);

/**
 * Update all graphics items to represent an updated conformer
 * @param atom_items All atom items to update
 * @param bond_items All bond items to update
 * @param sgroup_items All sgroup items to update
 * @param mol The molecule containing the updated conformer
 */
SKETCHER_API void
update_conf_for_mol_graphics_items(const QList<QGraphicsItem*>& atom_items,
                                   const QList<QGraphicsItem*>& bond_items,
                                   const QList<QGraphicsItem*>& sgroup_items,
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
 * @param arrow_color the color to use for the arrow itself
 * @param outline_color the color to use for the thin outline around the arrow
 * (used to ensure that the cursor doesn't disappear when its over something
 * that's the same color as the arrow)
 */
SKETCHER_API QPixmap get_arrow_cursor_pixmap(const QColor& arrow_color,
                                             const QColor& outline_color);

/**
 * Get a predictive hightlighting path for all atoms and bonds contained within
 * the given S-group
 * @param s_group The S-group to generate the path for
 * @return The newly generated highlighting path
 */
SKETCHER_API QPainterPath
get_predictive_highlighting_path_for_s_group_atoms_and_bonds(
    const RDKit::SubstanceGroup& s_group);

SKETCHER_API QPainterPath
get_selection_highlighting_path_for_atom(const RDKit::Atom* atom);

SKETCHER_API QPainterPath
get_predictive_highlighting_path_for_atom(const RDKit::Atom* atom);

SKETCHER_API QPainterPath
get_selection_highlighting_path_for_bond(const RDKit::Bond* bond);

SKETCHER_API QPainterPath
get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond);

/**
 * Calculate a painter path around the perimeter of `line` if `line` were drawn
 * with a pen width of `2 * half_width`.
 *
 * @param line The line for the path to go around
 * @param half_width The desired half-width for the painter path
 */
SKETCHER_API QPainterPath path_around_line(const QLineF& line,
                                           const qreal half_width);
/**
 * @return whether a graphics item is one of the types specified by the type
 * flag.
 */
SKETCHER_API bool item_matches_type_flag(QGraphicsItem* item,
                                         InteractiveItemFlagType type_flag);

/**
 * Return all model objects that are represented by the given collection of
 * graphics items.
 */
template <typename T>
SKETCHER_API std::tuple<std::unordered_set<const RDKit::Atom*>,
                        std::unordered_set<const RDKit::Bond*>,
                        std::unordered_set<const RDKit::SubstanceGroup*>,
                        std::unordered_set<const NonMolecularObject*>>
get_model_objects_for_graphics_items(const T& items);

} // namespace sketcher
} // namespace schrodinger
