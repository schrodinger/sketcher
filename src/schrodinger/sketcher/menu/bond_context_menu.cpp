#include "schrodinger/sketcher/menu/bond_context_menu.h"

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

ModifyBondsMenu::ModifyBondsMenu(QWidget* parent) : QMenu(parent)
{
    auto addBondAction = [this](auto text, auto type) {
        addAction(text, this,
                  [this, type]() { emit changeTypeRequested(type); });
    };

    setTitle("Modify Bonds");
    m_flip_action =
        addAction("Flip Substituent", this, &ModifyBondsMenu::flipRequested);
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
}

void ModifyBondsMenu::setFlipEnabled(bool b)
{
    m_flip_action->setEnabled(b);
}

QMenu* ModifyBondsMenu::createOtherTypeMenu()
{
    auto other_type_menu = new QMenu("Other Type", this);

    auto addBondAction = [this, other_type_menu](auto text, auto type) {
        other_type_menu->addAction(
            text, this, [this, type]() { emit changeTypeRequested(type); });
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
        query_menu->addAction(
            text, this, [this, type]() { emit changeTypeRequested(type); });
    };

    addBondAction("Any", BondTool::ANY);
    addBondAction("Single/Double", BondTool::SINGLE_OR_DOUBLE);
    addBondAction("Double/Aromatic", BondTool::DOUBLE_OR_AROMATIC);
    addBondAction("Single/Aromatic", BondTool::SINGLE_OR_AROMATIC);

    return query_menu;
}

BondContextMenu::BondContextMenu(QWidget* parent) : ModifyBondsMenu(parent)
{
    addSeparator();
    addAction("Delete", this, &BondContextMenu::deleteRequested);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/bond_context_menu.moc"
