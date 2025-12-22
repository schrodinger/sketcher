#include "schrodinger/sketcher/menu/bond_context_menu.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"

#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

ModifyBondsMenu::ModifyBondsMenu(QWidget* parent) : AbstractContextMenu(parent)
{
    auto addBondAction = [this](auto text, auto type) {
        addAction(text, this,
                  [this, type]() { emit changeTypeRequested(type, m_bonds); });
    };
    setTitle("Modify Bonds");
    m_flip_action = addAction("Flip Substituent", this,
                              [this]() { emit flipRequested(m_bonds); });

    addSeparator();
    addBondAction("Single", BondTool::SINGLE);
    addBondAction("Double", BondTool::DOUBLE);
    addBondAction("Triple", BondTool::TRIPLE);
    addBondAction("Aromatic", BondTool::AROMATIC);
    addSeparator();
    addBondAction("Up", BondTool::SINGLE_UP);
    addBondAction("Down", BondTool::SINGLE_DOWN);
    addMenu(createOtherTypeMenu());
    addMenu(createQueryMenu());
    addMenu(createTopologyMenu());
}

void ModifyBondsMenu::updateActions()
{
    // exactly 1 bonds should be selected and only non-ring bonds can be
    // flipped to avoid distorting the structure
    auto is_in_ring = [](auto bond) {
        auto& mol = bond->getOwningMol();
        return mol.getRingInfo()->numBondRings(bond->getIdx()) > 0;
    };
    auto flip_enabled = m_bonds.size() == 1 && !is_in_ring(*m_bonds.begin());
    m_flip_action->setEnabled(flip_enabled);
}

void ModifyBondsMenu::setFlipVisible(bool b)
{
    m_flip_action->setVisible(b);
}

QMenu* ModifyBondsMenu::createOtherTypeMenu()
{
    auto other_type_menu = new QMenu("Other Type", this);

    auto addBondAction = [this, other_type_menu](auto text, auto type) {
        other_type_menu->addAction(text, this, [this, type]() {
            emit changeTypeRequested(type, m_bonds);
        });
    };

    addBondAction("Coordinate", BondTool::COORDINATE);
    addBondAction("Zero Order", BondTool::ZERO);
    other_type_menu->addSeparator();
    addBondAction("Single Up/Down", BondTool::SINGLE_EITHER);
    addBondAction("Double Cis/Trans", BondTool::DOUBLE_EITHER);

    return other_type_menu;
}

QMenu* ModifyBondsMenu::createQueryMenu()
{
    auto query_menu = new QMenu("Query", this);

    auto addBondAction = [this, query_menu](auto text, auto type) {
        query_menu->addAction(text, this, [this, type]() {
            emit changeTypeRequested(type, m_bonds);
        });
    };

    addBondAction("Any", BondTool::ANY);
    addBondAction("Single/Double", BondTool::SINGLE_OR_DOUBLE);
    addBondAction("Double/Aromatic", BondTool::DOUBLE_OR_AROMATIC);
    addBondAction("Single/Aromatic", BondTool::SINGLE_OR_AROMATIC);

    return query_menu;
}

QMenu* ModifyBondsMenu::createTopologyMenu()
{
    auto topology_menu = new QMenu("Topology", this);

    auto addBondQueryAction = [this, topology_menu](auto text, auto type) {
        topology_menu->addAction(text, this, [this, type]() {
            emit changeQueryRequested(type, m_bonds);
        });
    };
    addBondQueryAction("In Ring", schrodinger::sketcher::BondTopology::IN_RING);
    addBondQueryAction("Not In a Ring",
                       schrodinger::sketcher::BondTopology::NOT_IN_RING);
    addBondQueryAction("Either", schrodinger::sketcher::BondTopology::EITHER);
    return topology_menu;
}

BondContextMenu::BondContextMenu(QWidget* parent) : ModifyBondsMenu(parent)
{
    addSeparator();
    addAction("Delete", this, [this]() {
        emit deleteRequested(m_bonds, m_secondary_connections);
    });
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/bond_context_menu.moc"
