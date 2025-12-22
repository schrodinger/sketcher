#include "schrodinger/sketcher/menu/selection_context_menu.h"

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/sgroup.h"
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
    m_clean_up_region_action = addAction(
        "Clean Up Region", this, &SelectionContextMenu::cleanUpRegionRequested);
    addSeparator();
    addAction(m_cut_copy_actions->m_cut_action);
    addAction(m_cut_copy_actions->m_copy_action);
    addMenu(m_cut_copy_actions->m_copy_as_menu);
    addSeparator();
    m_flip_action =
        addAction("Flip", this, &SelectionContextMenu::flipRequested);
    m_flip_molecule_menu = new QMenu("Flip Molecule");
    m_flip_molecule_menu->addAction(
        "Horizontally", this, &SelectionContextMenu::flipHorizontalRequested);
    m_flip_molecule_menu->addAction(
        "Vertically", this, &SelectionContextMenu::flipVerticalRequested);
    addMenu(m_flip_molecule_menu);
    addMenu(m_modify_atoms_menu);
    addMenu(m_modify_bonds_menu);
    addMenu(createAddToSelectionMenu());
    addSeparator();
    addAction("Delete", this, [this]() {
        emit deleteRequested(m_atoms, m_bonds, m_secondary_connections,
                             m_sgroups, m_non_molecular_objects);
    });

    connect(m_cut_copy_actions, &CutCopyActionManager::cutRequested, this,
            &SelectionContextMenu::cutRequested);
    connect(m_cut_copy_actions, &CutCopyActionManager::copyRequested, this,
            &SelectionContextMenu::copyRequested);
    connect(m_cut_copy_actions, &CutCopyActionManager::copyAsImageRequested,
            this, &SelectionContextMenu::copyAsImageRequested);
}

void SelectionContextMenu::setContextItems(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::Bond*>& secondary_connections,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    m_modify_atoms_menu->setContextItems(atoms, bonds, secondary_connections,
                                         sgroups, non_molecular_objects);
    m_modify_bonds_menu->setContextItems(atoms, bonds, secondary_connections,
                                         sgroups, non_molecular_objects);
    AbstractContextMenu::setContextItems(atoms, bonds, secondary_connections,
                                         sgroups, non_molecular_objects);
}

void SelectionContextMenu::updateActions()
{
    // enable clean up region only if there's one continuous region
    m_clean_up_region_action->setEnabled(
        is_contiguous_region(m_atoms, m_bonds));

    bool enable = m_atoms.size() >= 2 && in_same_fragment(m_atoms);
    m_variable_bond_action->setEnabled(enable);

    auto* mol = m_mol_model->getMol();
    bool enable_bracket = !m_atoms.empty() && m_sgroups.empty() &&
                          can_atoms_form_sgroup(m_atoms, *mol) &&
                          !get_existing_sgroup_for_atoms(m_atoms, *mol);
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
    bool one_fragment = !m_atoms.empty() && in_same_fragment(m_atoms);
    m_flip_action->setEnabled(true);
    m_flip_action->setVisible(crossing_bonds_count == 1 && one_fragment);
    m_flip_molecule_menu->menuAction()->setVisible(crossing_bonds_count == 0 &&
                                                   one_fragment);

    // if neither flip actions is displayed, show a disabled flip action
    if (!m_flip_action->isVisible() &&
        !m_flip_molecule_menu->menuAction()->isVisible()) {
        m_flip_action->setVisible(true);
        m_flip_action->setEnabled(false);
    }
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

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/selection_context_menu.moc"
