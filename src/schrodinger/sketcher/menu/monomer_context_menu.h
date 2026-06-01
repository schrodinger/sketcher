#pragma once

#include <unordered_set>
#include <vector>

#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"
#include "schrodinger/sketcher/rdkit/monomer_analog.h"

class QAction;
class QMenu;

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Context menu for monomer beads. Dispatches by monomer type:
 *   - PEPTIDE: Mutate Residue, Set D/L-Form, Protonate (stub), Delete
 *   - NA_BASE: Mutate Base (A/C/G/U/T + DB analogs), Delete
 *   - NA_SUGAR: Change R↔dR (direction from m_primary_atom), Delete
 *   - NA_PHOSPHATE / CHEM: Delete only
 */
class SKETCHER_API MonomerContextMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    MonomerContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested(const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * @brief Mutate one or more groups of monomers in a single user action.
     * @param mutations is a list of (atoms, target HELM symbol) pairs
     * @param description is the undo-stack label the receiver should use
     * when coalescing a multi-pair batch into one undo entry.
     */
    void mutateMonomerRequested(std::vector<MonomerMutation> mutations,
                                QString description);

  protected:
    void updateActions() override;

  private:
    void createMutateResidueSubMenu();
    void createSetDFormAction();
    void createProtonateAction();
    void createMutateBaseSubMenu();
    void createSugarToggleAction();
    void createDeleteAction();

    QMenu* m_mutate_residue_menu = nullptr;
    QAction* m_set_d_form_action = nullptr;
    QAction* m_protonate_action = nullptr;
    QMenu* m_mutate_base_menu = nullptr;
    QAction* m_sugar_toggle_action = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
