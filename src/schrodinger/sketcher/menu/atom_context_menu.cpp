#include "schrodinger/sketcher/menu/atom_context_menu.h"

#include <QWidgetAction>

#include "schrodinger/sketcher/Atom.h"
#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/atom_utils.h"
#include "schrodinger/sketcher/constants.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/widget/set_atom_widget.h"

namespace schrodinger
{
namespace sketcher
{

ModifyAtomsMenu::ModifyAtomsMenu(SketcherModel* model, QWidget* parent) :
    QMenu(parent),
    m_sketcher_model(model)
{
    // Create new sketcher model for the set atom menu widget.
    m_set_atom_model = new SketcherModel(this);
    // Assign non-atom draw tool so that none of the buttons in the set atoms
    // menu widget will be highlighted.
    m_set_atom_model->setValue(ModelKey::DRAW_TOOL, DrawTool::BOND);
    m_replace_with_menu = new ReplaceAtomsWithMenu(model, this);
    connect(m_replace_with_menu,
            &ReplaceAtomsWithMenu::
                showEditAtomPropertiesDialogWithAllowedListRequested,
            this, [this]() { emit showEditAtomPropertiesRequested(true); });
    connect(m_replace_with_menu, &ReplaceAtomsWithMenu::newRGroupRequested,
            this, &ModifyAtomsMenu::newRGroupRequested);
    connect(m_replace_with_menu, &ReplaceAtomsWithMenu::existingRGroupRequested,
            this, &ModifyAtomsMenu::existingRGroupRequested);

    setTitle("Modify Atoms");

    addMenu(createElementMenu());
    m_increase_charge_act =
        addAction("+ Charge", this, &ModifyAtomsMenu::increaseChargeRequested);
    m_decrease_charge_act =
        addAction("– Charge", this, &ModifyAtomsMenu::decreaseChargeRequested);
    addSeparator();
    m_add_remove_explicit_h_act =
        addAction("Add Explicit Hydrogens", this,
                  &ModifyAtomsMenu::addRemoveExplicitHydrogensRequested);
    m_add_unpaired_electrons_act =
        addAction("Add Unpaired Electron", this,
                  &ModifyAtomsMenu::addUnpairedElectronRequested);
    m_remove_unpaired_electrons_act =
        addAction("Remove Unpaired Electron", this,
                  &ModifyAtomsMenu::removeUnpairedElectronRequested);
    addSeparator();
    m_edit_atom_properties_act =
        addAction("Edit Atom Properties...", this,
                  [this]() { emit showEditAtomPropertiesRequested(false); });

    addMenu(m_replace_with_menu);

    connect(m_sketcher_model, &SketcherModel::selectionChanged, this,
            &ModifyAtomsMenu::updateActionsEnabled);
    connect(m_element_widget, &SetAtomMenuWidget::anyButtonClicked, this,
            &QMenu::close);
    connect(m_set_atom_model, &SketcherModel::valuePinged, this,
            &ModifyAtomsMenu::onSetAtomModelPinged);
}

void ModifyAtomsMenu::updateActionsEnabled()
{
    std::unordered_set<sketcherAtom*> atoms;
    std::vector<sketcherAtom*> element_atoms;
    std::vector<sketcherAtom*> atoms_that_can_have_charge;

    for (const auto& obj : m_sketcher_model->getContextMenuObjects()) {
        auto atom = dynamic_cast<sketcherAtom*>(obj);
        if (atom != nullptr) {
            atoms.insert(atom);
            if (is_atomic_number(atom->getAtomType())) {
                element_atoms.push_back(atom);
            }
            if (can_have_charge(*atom)) {
                atoms_that_can_have_charge.push_back(atom);
            }
        }
    }

    // Disable hydrogen and unpaired e for non elements
    if (element_atoms.empty()) {
        m_add_remove_explicit_h_act->setEnabled(false);
        m_remove_unpaired_electrons_act->setEnabled(false);
        m_add_unpaired_electrons_act->setEnabled(false);
    } else {
        // Toggle add/remove explicit hydrogens
        m_add_remove_explicit_h_act->setEnabled(true);
        if (std::any_of(element_atoms.begin(), element_atoms.end(),
                        has_implicit_hydrogens)) {
            m_add_remove_explicit_h_act->setText("Add Explicit Hydrogens");
        } else {
            m_add_remove_explicit_h_act->setText("Remove Explicit Hydrogens");
        }

        // Cap unpaired electrons
        std::vector<int> unpaired_e_counts;
        std::transform(element_atoms.begin(), element_atoms.end(),
                       std::back_inserter(unpaired_e_counts),
                       [](auto a) { return a->getUnpairedElectronsN(); });
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
                       [](auto a) { return a->getCharge(); });

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
        emit requestContextMenuElementChange(element);
    }
}

void ModifyAtomsMenu::showEvent(QShowEvent* event)
{
    updateActionsEnabled();
    QMenu::showEvent(event);
}

QMenu* ModifyAtomsMenu::createElementMenu()
{
    auto element_menu = new QMenu("Set Element", this);
    auto action = new QWidgetAction(element_menu);
    m_element_widget = new SetAtomMenuWidget(element_menu);
    m_element_widget->setModel(m_set_atom_model);
    action->setDefaultWidget(m_element_widget);
    element_menu->addAction(action);
    return element_menu;
}

ReplaceAtomsWithMenu::ReplaceAtomsWithMenu(SketcherModel* model,
                                           QWidget* parent) :
    QMenu(parent),
    m_sketcher_model(model)
{
    setTitle("Replace Atoms with");

    addMenu(createWildcardMenu());
    m_allowed_list_act =
        addAction("Allowed List...", this,
                  &ReplaceAtomsWithMenu::
                      showEditAtomPropertiesDialogWithAllowedListRequested);
    addSeparator();

    addAction("New R-Group", this, &ReplaceAtomsWithMenu::newRGroupRequested);
    m_existing_rgroup_menu = new ExistingRGroupMenu(model, this);
    addMenu(m_existing_rgroup_menu);

    connect(m_existing_rgroup_menu,
            &ExistingRGroupMenu::existingRGroupRequested, this,
            &ReplaceAtomsWithMenu::existingRGroupRequested);
}

void ReplaceAtomsWithMenu::showEvent(QShowEvent* event)
{
    auto rgroup_numbers = m_sketcher_model->getRGroupNumbers();
    m_existing_rgroup_menu->setEnabled(!rgroup_numbers.empty());
    QMenu::showEvent(event);
}

QMenu* ReplaceAtomsWithMenu::createWildcardMenu()
{
    auto wildcard_menu = new QMenu("Wildcard", this);

    auto addBondAction = [this, wildcard_menu](auto text, auto type) {
        wildcard_menu->addAction(
            text, this, [this, type]() { emit changeTypeRequested(type); });
    };

    addBondAction("A (Any heavy atom)", A_QUERY_KEY);
    addBondAction("Q (Heteroatom)", Q_QUERY_KEY);
    addBondAction("M (Metal)", M_QUERY_KEY);
    addBondAction("X (Halogen)", X_QUERY_KEY);
    wildcard_menu->addSeparator();
    addBondAction("AH (Any or H)", AH_QUERY_KEY);
    addBondAction("QH (Hetero or H)", QH_QUERY_KEY);
    addBondAction("MH (Metal or H)", MH_QUERY_KEY);
    addBondAction("XH (Halogen or H)", XH_QUERY_KEY);
    return wildcard_menu;
}

ExistingRGroupMenu::ExistingRGroupMenu(SketcherModel* model, QWidget* parent) :
    QMenu(parent),
    m_sketcher_model(model)
{
    setTitle("Existing R-Group");
}

void ExistingRGroupMenu::showEvent(QShowEvent* event)
{
    updateItems();
    QMenu::showEvent(event);
}

void ExistingRGroupMenu::updateItems()
{
    clear();
    auto rgroup_numbers = m_sketcher_model->getRGroupNumbers();
    for (auto rgroup_number : rgroup_numbers) {
        auto title = QString("R%1").arg(rgroup_number);
        addAction(title, [this, rgroup_number]() {
            emit existingRGroupRequested(rgroup_number);
        });
    }
}

AtomContextMenu::AtomContextMenu(SketcherModel* model, QWidget* parent) :
    ModifyAtomsMenu(model, parent)
{
    // Prepend 'Invert Selection' and separator
    auto top_action = actions().first();
    auto invert_selection_action = new QAction("Invert Selection", this);
    connect(invert_selection_action, &QAction::triggered, this,
            &AtomContextMenu::atomSelectionInvertRequested);
    insertAction(top_action, invert_selection_action);
    insertAction(top_action, addSeparator());

    // Rename replace menu to 'Replace with'
    m_replace_with_menu->setTitle("Replace with");

    m_add_brackets_act = new QAction("Add Brackets...", this);
    connect(m_add_brackets_act, &QAction::triggered, this,
            &AtomContextMenu::bracketSubgroupDialogRequested);
    insertAction(m_replace_with_menu->menuAction(), m_add_brackets_act);

    // Append 'Delete' action
    addSeparator();
    addAction("Delete", this, &AtomContextMenu::deleteRequested);

    // Show 'Invert Selection' only if there is an active selection
    invert_selection_action->setVisible(model->hasActiveSelection());
}

void AtomContextMenu::setAddToBracketGroupEnabled(bool enabled)
{
    m_add_brackets_act->setEnabled(enabled);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/atom_context_menu.moc"
