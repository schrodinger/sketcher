
#include "schrodinger/sketcher/menu/bracket_subgroup_context_menu.h"

namespace schrodinger
{
namespace sketcher
{

BracketSubgroupContextMenu::BracketSubgroupContextMenu(QWidget* parent) :
    AbstractContextMenu(parent)
{
    m_modify_notation_action = addAction("Modify Notation...", this, [this]() {
        emit bracketSubgroupDialogRequested(m_sgroups);
    });

    addAction("Remove Brackets", this,
              [this]() { emit deleteRequested(m_sgroups); });
}

void BracketSubgroupContextMenu::updateActions()
{
    m_modify_notation_action->setEnabled(m_sgroups.size() == 1);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/bracket_subgroup_context_menu.moc"
