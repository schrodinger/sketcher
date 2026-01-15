#pragma once

#include <unordered_set>

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
class Bond;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class NonMolecularObject;

class SKETCHER_API AbstractContextMenu : public QMenu
{
    Q_OBJECT
  public:
    AbstractContextMenu(QWidget* parent = nullptr);

    /**
     * Update the items associated with the context menu actions
     *
     * These items are pointers to components of the model as derived from
     * the graphical objects in the scene. We expect these context menus to:
     * 1. be populated with these references when shown
     * 2. block all other interactions with the scene while shown
     * 3. immediately hidden when an action is clicked
     * In other words, we expect that the model will not change for the
     * duration of the context menu being shown.
     */
    virtual void setContextItems(
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::Bond*>& secondary_connections,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::unordered_set<const NonMolecularObject*>&
            non_molecular_objects);

  protected:
    /**
     * Overrides the show event to also update enabled actions as appropriate
     */
    void showEvent(QShowEvent* event) override;

    /**
     * Pure virtual method required to update enabled menu actions
     */
    virtual void updateActions() = 0;

    std::unordered_set<const RDKit::Atom*> m_atoms;
    std::unordered_set<const RDKit::Bond*> m_bonds;
    std::unordered_set<const RDKit::Bond*> m_secondary_connections;
    std::unordered_set<const RDKit::SubstanceGroup*> m_sgroups;
    std::unordered_set<const NonMolecularObject*> m_non_molecular_objects;
};

} // namespace sketcher
} // namespace schrodinger