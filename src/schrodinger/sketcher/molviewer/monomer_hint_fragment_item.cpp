#include "schrodinger/sketcher/molviewer/monomer_hint_fragment_item.h"

#include <boost/range/join.hpp>

#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"

namespace schrodinger::sketcher
{

MonomerHintFragmentItem::MonomerHintFragmentItem(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const int atom_index_to_hide, const int bond_index_to_label,
    const QColor monomer_background_color, QGraphicsItem* parent) :
    QGraphicsItemGroup(parent),
    m_frag(fragment),
    m_fonts(&fonts),
    m_atom_index_to_hide(atom_index_to_hide),
    m_bond_index_to_label(bond_index_to_label),
    m_monomer_background_color(monomer_background_color)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER_FRAGMENT_HINT));
    createGraphicsItems();
}

void MonomerHintFragmentItem::createGraphicsItems()
{
    auto [all_items, atom_to_atom_item, bond_to_bond_item,
          bond_to_secondary_connection_item, s_group_to_s_group_item] =
        create_graphics_items_for_mol(&m_frag, *m_fonts);
    for (auto* item : all_items) {
        addToGroup(item);
    }
    for (auto& kv : atom_to_atom_item) {
        if (auto* monomer_item =
                dynamic_cast<AbstractMonomerItem*>(kv.second)) {
            // if we used a transparent monomer background color, then we'd be
            // able to see the bond behind the letter.  To avoid that, we use
            // the same color as the Scene's background.
            monomer_item->setMonomerColors(m_monomer_background_color,
                                           CURSOR_HINT_COLOR,
                                           CURSOR_HINT_COLOR);
        }
        // hide the monomer where this fragment connects to the existing
        // structure
        if (kv.first->getIdx() == m_atom_index_to_hide) {
            kv.second->setVisible(false);
        }
        m_atom_items.append(kv.second);
    }
    for (auto& kv : boost::range::join(bond_to_bond_item,
                                       bond_to_secondary_connection_item)) {
        if (auto* connector_item =
                qgraphicsitem_cast<MonomerConnectorItem*>(kv.second)) {
            connector_item->setConnectorStyle(
                CURSOR_HINT_COLOR, MONOMER_FRAGMENT_HINT_CONNECTOR_WIDTH);
        }
    }

    // label the attachment points
    if (m_bond_index_to_label >= 0) {
        auto* bond = m_frag.getBondWithIdx(m_bond_index_to_label);
        // TODO: If frag contains connections with arrowheads, this function
        //       needs the monomer items so it can figure out the arrowhead
        //       offset. These would normally be retrieved using the Scene, but
        //       that won't work for hint structures since they're not included
        //       in the Scene's bookkeeping. For now, we just pass nullptr for
        //       Scene and avoid fragments with arrowheads.
        auto items = create_attachment_point_labels_for_connector(
            bond, false, STRUCTURE_HINT_COLOR, *m_fonts, nullptr);
        for (auto* item : items) {
            addToGroup(item);
        }
    }
}

} // namespace schrodinger::sketcher
