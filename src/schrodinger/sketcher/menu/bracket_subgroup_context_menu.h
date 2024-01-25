#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

namespace RDKit
{
class SubstanceGroup;
} // namespace RDKit
namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API BracketSubgroupContextMenu : public AbstractContextMenu
{
    Q_OBJECT

  public:
    BracketSubgroupContextMenu(QWidget* parent = nullptr);

  signals:
    void bracketSubgroupDialogRequested(
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups) const;
    void deleteRequested(
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups) const;

  protected:
    virtual void updateActions();

    QAction* m_modify_notation_action = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
