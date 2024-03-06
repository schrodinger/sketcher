#include "schrodinger/sketcher/menu/selection_context_menu.h"

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/sgroup.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/menu/bond_context_menu.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/subset.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

SelectionContextMenu::SelectionContextMenu(SketcherModel* model,
                                           MolModel* mol_model,
                                           QWidget* parent) :
    AbstractContextMenu(parent),
    m_sketcher_model(model),
    m_mol_model(mol_model)
{
    m_cut_copy_actions = new CutCopyActionManager(this);
    m_cut_copy_actions->setModel(model);
    m_modify_atoms_menu = new ModifyAtomsMenu(model, mol_model, this);
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
    addMenu(createAddToSelectionMenu());
    addMenu(createReplaceSelectionWithMenu());
    addSeparator();
    addAction("Delete", this, [this]() {
        emit deleteRequested(m_atoms, m_bonds, m_sgroups,
                             m_non_molecular_objects);
    });

    connect(m_cut_copy_actions, &CutCopyActionManager::cutRequested, this,
            &SelectionContextMenu::cutRequested);
    connect(m_cut_copy_actions, &CutCopyActionManager::copyRequested, this,
            &SelectionContextMenu::copyRequested);
}

void SelectionContextMenu::setContextItems(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    m_modify_atoms_menu->setContextItems(atoms, bonds, sgroups,
                                         non_molecular_objects);
    m_modify_bonds_menu->setContextItems(atoms, bonds, sgroups,
                                         non_molecular_objects);
    AbstractContextMenu::setContextItems(atoms, bonds, sgroups,
                                         non_molecular_objects);
}

void SelectionContextMenu::updateActions()
{
    bool enable = m_atoms.size() >= 2 && in_same_fragment(m_atoms);
    m_variable_bond_action->setEnabled(enable);

    auto rgroup_numbers = get_all_r_group_numbers(m_mol_model->getMol());
    m_existing_rgroup_menu->setEnabled(!rgroup_numbers.empty());

    auto* mol = m_mol_model->getMol();
    bool enable_bracket =
        !m_atoms.empty() && m_sgroups.empty() &&
        rdkit_extensions::can_atoms_form_sgroup(m_atoms, *mol) &&
        !rdkit_extensions::get_existing_sgroup_for_atoms(m_atoms, *mol);
    m_bracket_group_action->setEnabled(enable_bracket);

    /**
     * enable flip action if the selection is in one fragment and there is
     * exactly one bond that connects selected atoms to atoms outside the
     * selection
     */
    auto all_bonds = mol->bonds();
    auto crossing_bonds_count =
        std::count_if(all_bonds.begin(), all_bonds.end(), [this](auto bond) {
            return m_atoms.count(bond->getBeginAtom()) !=
                   m_atoms.count(bond->getEndAtom());
        });
    m_flip_action->setEnabled(crossing_bonds_count == 1 &&
                              in_same_fragment(m_atoms));
}

QMenu* SelectionContextMenu::createAddToSelectionMenu()
{
    auto add_to_selection_menu = new QMenu("Add to Selection", this);
    m_bracket_group_action =
        add_to_selection_menu->addAction("Bracket Subgroup...", this, [this]() {
            emit bracketSubgroupDialogRequested(m_atoms);
        });
    m_variable_bond_action = add_to_selection_menu->addAction(
        "Variable Attachment Bond", this,
        [this]() { emit variableAttachmentBondRequested(m_atoms); });
    return add_to_selection_menu;
}

QMenu* SelectionContextMenu::createReplaceSelectionWithMenu()
{
    auto replace_selection_menu = new QMenu("Replace Selection with", this);
    replace_selection_menu->addAction(
        "New R-Group", this, [this]() { emit newRGroupRequested(m_atoms); });

    m_existing_rgroup_menu =
        new ExistingRGroupMenu(m_mol_model, replace_selection_menu);
    connect(m_existing_rgroup_menu,
            &ExistingRGroupMenu::existingRGroupRequested, this,
            [this](auto rgroup_number) {
                emit existingRGroupRequested(m_atoms, rgroup_number);
            });
    replace_selection_menu->addMenu(m_existing_rgroup_menu);

    // SKETCH-951 TODO: add QAction for "Aromatized Structure"
    // SKETCH-951 TODO: add QAction for "Kekulized Structure"
    return replace_selection_menu;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/selection_context_menu.moc"
