#include "schrodinger/sketcher/menu/selection_context_menu.h"

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/menu/bond_context_menu.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

SelectionContextMenu::SelectionContextMenu(SketcherModel* model,
                                           QWidget* parent) :
    QMenu(parent),
    m_sketcher_model(model)
{
    m_cut_copy_actions = new CutCopyActionManager(this);
    m_cut_copy_actions->setModel(model);
    m_modify_atoms_menu = new ModifyAtomsMenu(model, this);
    m_modify_bonds_menu = new ModifyBondsMenu(this);
    m_modify_bonds_menu->setFlipVisible(false);

    addAction("Invert Selection", this,
              &SelectionContextMenu::invertSelectionRequested);
    addSeparator();
    addAction(m_cut_copy_actions->m_cut_action);
    addAction(m_cut_copy_actions->m_copy_action);
    addMenu(m_cut_copy_actions->m_copy_as_menu);
    addSeparator();
    m_flip_action =
        addAction("Flip", this, &SelectionContextMenu::flipRequested);
    addMenu(m_modify_atoms_menu);
    addMenu(m_modify_bonds_menu);
    addMenu(createAddToSelectionMenu(model));
    addMenu(createReplaceSelectionWithMenu(model));
    addSeparator();
    addAction("Delete", this, &SelectionContextMenu::deleteRequested);

    connect(m_cut_copy_actions, &CutCopyActionManager::cutRequested, this,
            &SelectionContextMenu::cutRequested);
    connect(m_cut_copy_actions, &CutCopyActionManager::copyRequested, this,
            &SelectionContextMenu::copyRequested);
}

void SelectionContextMenu::showEvent(QShowEvent* event)
{
    QMenu::showEvent(event);
    updateActionsEnabled();
}

void SelectionContextMenu::updateActionsEnabled()
{
    auto same_mol = [](const auto& atoms) {
        auto mol = atoms.front()->getOwningMol();
        return std::all_of(atoms.begin(), atoms.end(), [mol](auto atom) {
            return &atom->getOwningMol() == &mol;
        });
    };

    auto atoms = m_sketcher_model->getContextMenuAtoms();
    bool enable = atoms.size() >= 2 && same_mol(atoms);
    m_variable_bond_action->setEnabled(enable);
}

void SelectionContextMenu::setConvertToBracketGroupEnabled(bool b)
{
    m_bracket_group_action->setEnabled(b);
}

void SelectionContextMenu::setFlipEnabled(bool b)
{
    m_flip_action->setEnabled(b);
}

QMenu* SelectionContextMenu::createAddToSelectionMenu(SketcherModel* model)
{
    auto add_to_selection_menu = new QMenu("Add to Selection", this);
    m_bracket_group_action = add_to_selection_menu->addAction(
        "Bracket Subgroup...", this,
        &SelectionContextMenu::bracketSubgroupDialogRequested);
    m_variable_bond_action = add_to_selection_menu->addAction(
        "Variable Attachment Bond", this,
        &SelectionContextMenu::variableAttachmentBondRequested);
    return add_to_selection_menu;
}

QMenu*
SelectionContextMenu::createReplaceSelectionWithMenu(SketcherModel* model)
{
    auto replace_selection_menu = new QMenu("Replace Selection with", this);

    auto existing_rgroup_menu =
        new ExistingRGroupMenu(model, replace_selection_menu);

    connect(existing_rgroup_menu, &ExistingRGroupMenu::existingRGroupRequested,
            this, &SelectionContextMenu::existingRGroupRequested);

    replace_selection_menu->addAction(
        "New R-Group", this, &SelectionContextMenu::newRGroupRequested);
    replace_selection_menu->addMenu(existing_rgroup_menu);
    // SKETCH-951 TODO: add QAction for "Aromatized Structure"
    // SKETCH-951 TODO: add QAction for "Kekulized Structure"
    return replace_selection_menu;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/selection_context_menu.moc"
