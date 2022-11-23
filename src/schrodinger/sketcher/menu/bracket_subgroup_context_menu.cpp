#include "schrodinger/sketcher/menu/bracket_subgroup_context_menu.h"

#include "schrodinger/rdkit_extensions/convert.h"

namespace schrodinger
{
namespace sketcher
{

BracketSubgroupContextMenu::BracketSubgroupContextMenu(QWidget* parent) :
    QMenu(parent)
{
    addAction("Modify Notation...", this,
              &BracketSubgroupContextMenu::bracketSubgroupDialogRequested);
    addAction("Remove Brackets", this,
              &BracketSubgroupContextMenu::deleteRequested);
}

} // namespace sketcher
} // namespace schrodinger

#include "bracket_subgroup_context_menu.moc"
