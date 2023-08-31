#include "schrodinger/sketcher/molviewer/scene_utils.h"

#include <QGraphicsItem>

#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

std::tuple<std::vector<QGraphicsItem*>,
           std::unordered_map<const RDKit::Atom*, AtomItem*>,
           std::unordered_map<const RDKit::Bond*, BondItem*>>
create_graphics_items_for_mol(const RDKit::ROMol* mol, const Fonts& fonts,
                              AtomItemSettings& atom_item_settings,
                              BondItemSettings& bond_item_settings,
                              const bool draw_attachment_points)
{
    unsigned int num_atoms = mol->getNumAtoms();
    if (num_atoms == 0) {
        // If there are no atoms, then there's nothing more to do.  Also, if
        // there are no atoms, then shorten_attachment_point_bonds() will raise
        // a ConformerException.
        return {{}, {}, {}};
    }
    RDKit::Conformer conformer = shorten_attachment_point_bonds(mol);

    std::vector<QGraphicsItem*> all_items;
    std::unordered_map<const RDKit::Atom*, AtomItem*> atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, BondItem*> bond_to_bond_item;

    // create atom items
    for (std::size_t i = 0; i < num_atoms; ++i) {
        const auto* atom = mol->getAtomWithIdx(i);
        if (!draw_attachment_points && is_attachment_point(atom)) {
            continue;
        }
        const auto pos = conformer.getAtomPos(i);
        auto* atom_item = new AtomItem(atom, fonts, atom_item_settings);
        atom_item->setPos(to_scene_xy(pos));
        atom_to_atom_item[atom] = atom_item;
        all_items.push_back(atom_item);
    }

    // create bond items
    for (auto bond : mol->bonds()) {
        if (!draw_attachment_points && is_attachment_point_bond(bond)) {
            continue;
        }
        const auto* from_atom_item = atom_to_atom_item[bond->getBeginAtom()];
        const auto* to_atom_item = atom_to_atom_item[bond->getEndAtom()];
        auto* bond_item = new BondItem(bond, *from_atom_item, *to_atom_item,
                                       fonts, bond_item_settings);
        bond_to_bond_item[bond] = bond_item;
        all_items.push_back(bond_item);
    }

    return {all_items, atom_to_atom_item, bond_to_bond_item};
}

void update_conf_for_mol_graphics_items(const QList<QGraphicsItem*>& atom_items,
                                        const QList<QGraphicsItem*>& bond_items,
                                        const RDKit::ROMol& mol)
{
    auto conf = mol.getConformer();
    for (auto* item : atom_items) {
        auto* atom_item = qgraphicsitem_cast<AtomItem*>(item);
        auto pos = conf.getAtomPos(atom_item->getAtom()->getIdx());
        atom_item->setPos(to_scene_xy(pos));
        atom_item->updateCachedData();
    }
    for (auto* item : bond_items) {
        auto* bond_item = qgraphicsitem_cast<BondItem*>(item);
        bond_item->updateCachedData();
    }
}

} // namespace sketcher
} // namespace schrodinger
