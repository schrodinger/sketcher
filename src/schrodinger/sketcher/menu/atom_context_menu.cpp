#include "schrodinger/sketcher/menu/atom_context_menu.h"

#include <QWidgetAction>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/constants.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/widget/set_atom_widget.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
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

void ModifyAtomsMenu::setContextAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    m_atoms = atoms;
    updateActionsEnabled();
}

void ModifyAtomsMenu::showEvent(QShowEvent* event)
{
    updateActionsEnabled();
    QMenu::showEvent(event);
}

void ModifyAtomsMenu::updateActionsEnabled()
{
    auto atoms = m_atoms;

    std::vector<const RDKit::Atom*> element_atoms;
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(element_atoms),
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
        if (std::any_of(element_atoms.begin(), element_atoms.end(),
                        [](auto a) { return a->getTotalNumHs() > 0; })) {
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
    // Rename replace menu to 'Replace with'
    m_replace_with_menu->setTitle("Replace with");

    m_add_brackets_act = new QAction("Add Brackets...", this);
    connect(m_add_brackets_act, &QAction::triggered, this,
            &AtomContextMenu::bracketSubgroupDialogRequested);
    insertAction(m_replace_with_menu->menuAction(), m_add_brackets_act);

    // Append 'Delete' action
    addSeparator();
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

void AtomContextMenu::updateActionsEnabled()
{
    auto has_two_bonds = [](auto atom) {
        auto [begin, end] = atom->getOwningMol().getAtomBonds(atom);
        return std::distance(begin, end) == 2;
    };
    bool enable = m_atoms.size() == 1 && has_two_bonds(*m_atoms.begin());
    m_add_brackets_act->setEnabled(enable);

    ModifyAtomsMenu::updateActionsEnabled();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/atom_context_menu.moc"
