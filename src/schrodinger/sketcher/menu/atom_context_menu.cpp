#include "schrodinger/sketcher/menu/atom_context_menu.h"

#include <QWidgetAction>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/rdkit/sgroup.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/widget/set_atom_widget.h"

namespace schrodinger
{
namespace sketcher
{

ModifyAtomsMenu::ModifyAtomsMenu(SketcherModel* model, MolModel* mol_model,
                                 QWidget* parent) :
    AbstractContextMenu(parent),
    m_sketcher_model(model)
{
    // Create new sketcher model for the set atom menu widget.
    m_set_atom_model = new SketcherModel(this);
    // Assign non-atom draw tool so that none of the buttons in the set atoms
    // menu widget will be highlighted.
    m_set_atom_model->setValue(ModelKey::DRAW_TOOL, DrawTool::BOND);
    m_replace_with_menu = new ReplaceAtomsWithMenu(mol_model, this);
    connect(m_replace_with_menu,
            &ReplaceAtomsWithMenu::
                showEditAtomPropertiesDialogWithAllowedListRequested,
            this, [this]() {
                emit showEditAtomPropertiesRequested(*m_atoms.begin(), true);
            });
    connect(m_replace_with_menu, &ReplaceAtomsWithMenu::newRGroupRequested,
            this, [this]() { emit newRGroupRequested(m_atoms); });
    connect(m_replace_with_menu, &ReplaceAtomsWithMenu::existingRGroupRequested,
            this, [this](auto rgroup_number) {
                emit existingRGroupRequested(m_atoms, rgroup_number);
            });
    connect(m_replace_with_menu, &ReplaceAtomsWithMenu::changeTypeRequested,
            this,
            [this](auto type) { emit changeTypeRequested(m_atoms, type); });
    setTitle("Modify Atoms");

    addMenu(createElementMenu());
    m_increase_charge_act = addAction("+ Charge", this, [this]() {
        emit adjustChargeRequested(m_atoms, +1);
    });
    m_decrease_charge_act = addAction("– Charge", this, [this]() {
        emit adjustChargeRequested(m_atoms, -1);
    });
    addSeparator();
    m_add_remove_explicit_h_act =
        addAction("Add Explicit Hydrogens", this, [this]() {
            emit addRemoveExplicitHydrogensRequested(m_atoms);
        });
    m_add_unpaired_electrons_act =
        addAction("Add Unpaired Electron", this, [this]() {
            emit adjustRadicalElectronsRequested(m_atoms, +1);
        });
    m_remove_unpaired_electrons_act =
        addAction("Remove Unpaired Electron", this, [this]() {
            emit adjustRadicalElectronsRequested(m_atoms, -1);
        });
    addSeparator();
    m_edit_atom_properties_act =
        addAction("Edit Atom Properties...", this, [this]() {
            emit showEditAtomPropertiesRequested(*m_atoms.begin(), false);
        });

    addMenu(m_replace_with_menu);
    connect(m_set_atom_model, &SketcherModel::valuePinged, this,
            &ModifyAtomsMenu::onSetAtomModelPinged);
}

void ModifyAtomsMenu::updateActions()
{
    auto atoms = m_atoms;

    std::unordered_set<const RDKit::Atom*> element_atoms;
    std::copy_if(atoms.begin(), atoms.end(),
                 std::inserter(element_atoms, element_atoms.begin()),
                 [](auto a) { return !a->hasQuery(); });

    std::vector<const RDKit::Atom*> atoms_that_can_have_charge;
    std::copy_if(atoms.begin(), atoms.end(),
                 std::back_inserter(atoms_that_can_have_charge),
                 [](auto a) { return !is_r_group(a); });

    // Disable hydrogen and unpaired e for non elements
    if (element_atoms.empty()) {
        m_add_remove_explicit_h_act->setEnabled(false);
        m_remove_unpaired_electrons_act->setEnabled(false);
        m_add_unpaired_electrons_act->setEnabled(false);
    } else {
        // Toggle add/remove explicit hydrogens
        m_add_remove_explicit_h_act->setEnabled(true);
        if (has_any_implicit_Hs(element_atoms)) {
            m_add_remove_explicit_h_act->setText("Add Explicit Hydrogens");
        } else {
            m_add_remove_explicit_h_act->setText("Remove Explicit Hydrogens");
        }

        // Cap unpaired electrons
        std::vector<int> unpaired_e_counts;
        std::transform(element_atoms.begin(), element_atoms.end(),
                       std::back_inserter(unpaired_e_counts),
                       [](auto a) { return a->getNumRadicalElectrons(); });
        const auto [min_unpaired_e_count, max_unpaired_e_count] =
            std::minmax_element(unpaired_e_counts.begin(),
                                unpaired_e_counts.end());
        m_remove_unpaired_electrons_act->setEnabled(*min_unpaired_e_count >
                                                    MIN_UNPAIRED_E);
        m_add_unpaired_electrons_act->setEnabled(*min_unpaired_e_count <
                                                 MAX_UNPAIRED_E);
    }

    if (atoms_that_can_have_charge.empty()) {
        // Disable if charges not allowed
        m_decrease_charge_act->setEnabled(false);
        m_increase_charge_act->setEnabled(false);
    } else {
        // Cap charges
        std::vector<int> charges;
        std::transform(atoms_that_can_have_charge.begin(),
                       atoms_that_can_have_charge.end(),
                       std::back_inserter(charges),
                       [](auto a) { return a->getFormalCharge(); });

        const auto [min_charge, max_charge] =
            std::minmax_element(charges.begin(), charges.end());

        m_decrease_charge_act->setEnabled(*min_charge > -ATOM_CHARGE_LIMIT);
        m_increase_charge_act->setEnabled(*max_charge < ATOM_CHARGE_LIMIT);
    }

    bool enable = atoms.size() == 1;
    m_edit_atom_properties_act->setEnabled(enable);
    m_replace_with_menu->m_allowed_list_act->setEnabled(enable);
}

void ModifyAtomsMenu::onSetAtomModelPinged(ModelKey key, QVariant value)
{
    if (key == ModelKey::ELEMENT) {
        auto element = Element(value.toInt());
        emit requestElementChange(m_atoms, element);
    }
}

QMenu* ModifyAtomsMenu::createElementMenu()
{
    auto element_menu = new QMenu("Set Element", this);
    auto action = new QWidgetAction(element_menu);
    auto element_widget = new SetAtomMenuWidget(element_menu);
    element_widget->setModel(m_set_atom_model);
    action->setDefaultWidget(element_widget);
    element_menu->addAction(action);
    connect(element_widget, &SetAtomMenuWidget::anyButtonClicked, action,
            &QAction::trigger);
    return element_menu;
}

ReplaceAtomsWithMenu::ReplaceAtomsWithMenu(MolModel* mol_model,
                                           QWidget* parent) :
    AbstractContextMenu(parent),
    m_mol_model(mol_model)
{
    setTitle("Replace Atoms with");

    addMenu(createWildcardMenu());
    m_allowed_list_act =
        addAction("Allowed List...", this,
                  &ReplaceAtomsWithMenu::
                      showEditAtomPropertiesDialogWithAllowedListRequested);
    addSeparator();

    addAction("New R-Group", this, &ReplaceAtomsWithMenu::newRGroupRequested);
    m_existing_rgroup_menu = new ExistingRGroupMenu(mol_model, this);
    addMenu(m_existing_rgroup_menu);

    connect(m_existing_rgroup_menu,
            &ExistingRGroupMenu::existingRGroupRequested, this,
            &ReplaceAtomsWithMenu::existingRGroupRequested);
}

void ReplaceAtomsWithMenu::updateActions()
{
    auto rgroup_numbers = get_all_r_group_numbers(m_mol_model->getMol());
    m_existing_rgroup_menu->setEnabled(!rgroup_numbers.empty());
}

QMenu* ReplaceAtomsWithMenu::createWildcardMenu()
{
    auto wildcard_menu = new QMenu("Wildcard", this);

    auto addQueryAction = [this, wildcard_menu](auto text, auto type) {
        wildcard_menu->addAction(
            text, this, [this, type]() { emit changeTypeRequested(type); });
    };

    addQueryAction("A (Any heavy atom)", AtomQuery::A);
    addQueryAction("Q (Heteroatom)", AtomQuery::Q);
    addQueryAction("M (Metal)", AtomQuery::M);
    addQueryAction("X (Halogen)", AtomQuery::X);
    wildcard_menu->addSeparator();
    addQueryAction("AH (Any or H)", AtomQuery::AH);
    addQueryAction("QH (Hetero or H)", AtomQuery::QH);
    addQueryAction("MH (Metal or H)", AtomQuery::MH);
    addQueryAction("XH (Halogen or H)", AtomQuery::XH);
    return wildcard_menu;
}

ExistingRGroupMenu::ExistingRGroupMenu(MolModel* mol_model, QWidget* parent) :
    AbstractContextMenu(parent),
    m_mol_model(mol_model)
{
    setTitle("Existing R-Group");
}

void ExistingRGroupMenu::updateActions()
{
    clear();
    auto rgroup_numbers = get_all_r_group_numbers(m_mol_model->getMol());
    for (auto rgroup_number : rgroup_numbers) {
        auto title = QString("R%1").arg(rgroup_number);
        addAction(title, [this, rgroup_number]() {
            emit existingRGroupRequested(rgroup_number);
        });
    }
}

AtomContextMenu::AtomContextMenu(SketcherModel* model, MolModel* mol_model,
                                 QWidget* parent) :
    ModifyAtomsMenu(model, mol_model, parent)
{
    // Rename replace menu to 'Replace with'
    m_replace_with_menu->setTitle("Replace with");

    m_add_brackets_act = new QAction("Add Brackets...", this);
    connect(m_add_brackets_act, &QAction::triggered, this,
            [this]() { emit bracketSubgroupDialogRequested(m_atoms); });
    insertAction(m_replace_with_menu->menuAction(), m_add_brackets_act);

    // Append 'Delete' action
    addSeparator();
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

void AtomContextMenu::updateActions()
{
    auto& mol = (*m_atoms.begin())->getOwningMol();
    bool enable_bracket = can_atoms_form_sgroup(m_atoms, mol) &&
                          !get_existing_sgroup_for_atoms(m_atoms, mol);
    m_add_brackets_act->setEnabled(enable_bracket);

    ModifyAtomsMenu::updateActions();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/atom_context_menu.moc"
