#pragma once

#include <unordered_set>

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

namespace RDKit
{
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

enum class BondTool;
enum class BondTopology;

class SKETCHER_API ModifyBondsMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    ModifyBondsMenu(QWidget* parent = nullptr);
    /**
     * Set whether the flip action is visible.
     */
    void setFlipVisible(bool b);

  signals:
    void
    changeTypeRequested(BondTool bond_tool,
                        const std::unordered_set<const RDKit::Bond*>& bonds);
    void
    changeQueryRequested(BondTopology bond_tool,
                         const std::unordered_set<const RDKit::Bond*>& bonds);
    void flipRequested(const std::unordered_set<const RDKit::Bond*>& bonds);

  protected:
    QAction* m_flip_action = nullptr;

  private:
    /**
     * Disable menu actions that don't make sense given the selected bonds
     */
    virtual void updateActions();

    QMenu* createOtherTypeMenu();
    QMenu* createQueryMenu();
    QMenu* createTopologyMenu();
};

class BondContextMenu : public ModifyBondsMenu
{
    Q_OBJECT
  public:
    BondContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested(
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::Bond*>& secondary_connections);
};

} // namespace sketcher
} // namespace schrodinger
