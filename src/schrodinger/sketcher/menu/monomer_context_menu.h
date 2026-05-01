#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

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
 * One unit of work for `MolModel::mutateMonomers`: a set of atoms that
 * should all be mutated to the same target HELM symbol. Used as the
 * batch element in `MonomerContextMenu::mutateResidueRequested`.
 */
struct MonomerMutation {
    std::unordered_set<const RDKit::Atom*> atoms;
    std::string helm_symbol;
};

/**
 * Context menu for monomer beads.
 *
 * For amino-acid (PEPTIDE) selections this exposes "Mutate Residue >"
 * with the natural AAs and their non-natural analogs, a dynamic
 * "Set D-Form / Set L-Form" toggle, a disabled "Protonate" stub, and
 * Delete. Other monomer types only see Delete.
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
    void mutateResidueRequested(std::vector<MonomerMutation> mutations,
                                QString description);

  protected:
    void updateActions() override;

  private:
    void createMutateResidueSubMenu();
    void createSetDFormAction();
    void createProtonateAction();
    void createDeleteAction();

    QMenu* m_mutate_residue_menu = nullptr;
    QAction* m_set_d_form_action = nullptr;
    QAction* m_protonate_action = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
