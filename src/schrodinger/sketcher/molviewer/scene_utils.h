#pragma once

#include <vector>
#include <unordered_map>

#include <QList>

#include "schrodinger/sketcher/definitions.h"

class QGraphicsItem;

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

} // namespace sketcher
} // namespace schrodinger
