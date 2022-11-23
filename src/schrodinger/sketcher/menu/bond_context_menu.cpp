#include "schrodinger/sketcher/menu/bond_context_menu.h"

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
    addBondAction("Single", sketcherBond::BondType::SINGLE);
    addBondAction("Double", sketcherBond::BondType::DOUBLE);
    addBondAction("Triple", sketcherBond::BondType::TRIPLE);
    addBondAction("Aromatic", sketcherBond::BondType::AROMATIC);
    addSeparator();
    addBondAction("Up", sketcherBond::BondType::SINGLE_UP);
    addBondAction("Down", sketcherBond::BondType::SINGLE_DOWN);
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

    addBondAction("Coordinate", sketcherBond::BondType::COORDINATE);
    addBondAction("Zero Order", sketcherBond::BondType::ZERO);
    other_type_menu->addSeparator();
    addBondAction("Single Up/Down", sketcherBond::BondType::SINGLE_EITHER);
    addBondAction("Double Cis/Trans", sketcherBond::BondType::DOUBLE_EITHER);

    return other_type_menu;
}

QMenu* ModifyBondsMenu::createQueryMenu()
{
    auto query_menu = new QMenu("Query", this);

    auto addBondAction = [this, query_menu](auto text, auto type) {
        query_menu->addAction(
            text, this, [this, type]() { emit changeTypeRequested(type); });
    };

    addBondAction("Any", sketcherBond::BondType::ANY);
    addBondAction("Single/Double", sketcherBond::BondType::SINGLE_OR_DOUBLE);
    addBondAction("Double/Aromatic",
                  sketcherBond::BondType::DOUBLE_OR_AROMATIC);
    addBondAction("Single/Aromatic",
                  sketcherBond::BondType::SINGLE_OR_AROMATIC);

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
