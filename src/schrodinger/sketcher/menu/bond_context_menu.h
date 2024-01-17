#pragma once

#include <unordered_set>

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

enum class BondTool;

class SKETCHER_API ModifyBondsMenu : public QMenu
{
    Q_OBJECT
  public:
    ModifyBondsMenu(QWidget* parent = nullptr);

    /**
     * Update the bonds associated with the context menu actions, updating
     * enabled actions as appropriate.
     *
     * These items are pointers to components of the model as derived from
     * the graphical objects in the scene. We expect these context menus to:
     * 1. be populated with these references when shown
     * 2. block all other interactions with the scene while shown
     * 3. immediately hidden when an action is clicked
     * In other words, we expect that the model will not change for the
     * duration of the context menu being shown.
     */
    void setContextBonds(const std::unordered_set<const RDKit::Bond*>& bonds);

    /**
     * Set whether the flip action is visible.
     */
    void setFlipVisible(bool b);
    // TODO: remove when enabling molviewer
    void setFlipEnabled(bool b);

  signals:
    void
    changeTypeRequested(BondTool bond_tool,
                        const std::unordered_set<const RDKit::Bond*>& bonds);
    void flipRequested(const std::unordered_set<const RDKit::Bond*>& bonds);

  protected:
    std::unordered_set<const RDKit::Bond*> m_bonds;

  private:
    QMenu* createOtherTypeMenu();
    QMenu* createQueryMenu();
    QAction* m_flip_action = nullptr;
};

class BondContextMenu : public ModifyBondsMenu
{
    Q_OBJECT
  public:
    BondContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested(const std::unordered_set<const RDKit::Bond*>& bonds);
};

} // namespace sketcher
} // namespace schrodinger
