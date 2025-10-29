#include "schrodinger/sketcher/widget/set_atom_widget.h"

#include <boost/assign.hpp>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_set_atom_widget.h"
#include "schrodinger/sketcher/widget/atom_query_popup.h"
#include "schrodinger/sketcher/widget/periodic_table_widget.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

using schrodinger::sketcher::AtomQuery;

namespace schrodinger
{
namespace sketcher
{

SetAtomWidget::SetAtomWidget(QWidget* parent) : AbstractDrawToolWidget(parent)
{
    ui.reset(new Ui::SetAtomWidget());
    ui->setupUi(this);

    m_atom_query_wdg = new AtomQueryPopup(this);
    ui->atom_query_btn->setEnumItem(static_cast<int>(AtomQuery::A));
    ui->atom_query_btn->setPopupWidget(m_atom_query_wdg);
    ui->last_picked_element_btn->setElement(Element::SI);

    for (auto btn : ui->atom_group->buttons()) {
        btn->setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);
    }
    ui->atom_query_btn->setStyleSheet(ATOM_QUERY_STYLE);

    m_periodic_table_wdg = new PeriodicTableWidget(this);
    m_periodic_table_wdg->setWindowFlags(Qt::Popup);
    ui->periodic_table_btn->setPopupDelay(0);
    ui->periodic_table_btn->setPopupWidget(m_periodic_table_wdg);
    ui->periodic_table_btn->showPopupIndicator(false);

    using ButtonElementBimapType = boost::bimap<QAbstractButton*, Element>;
    m_button_element_bimap =
        boost::assign::list_of<ButtonElementBimapType::relation>(
            ui->c_btn, Element::C)(ui->h_btn, Element::H)(
            ui->n_btn, Element::N)(ui->o_btn, Element::O)(
            ui->p_btn, Element::P)(ui->s_btn, Element::S)(
            ui->f_btn, Element::F)(ui->cl_btn, Element::CL);

    connect(ui->atom_group,
            static_cast<void (QButtonGroup::*)(QAbstractButton*)>(
                &QButtonGroup::buttonClicked),
            this, &SetAtomWidget::onAtomButtonClicked);
}

SetAtomWidget::~SetAtomWidget()
{
    delete m_atom_query_wdg;
}

void SetAtomWidget::setModel(SketcherModel* model)
{
    AbstractDrawToolWidget::setModel(model);
    m_periodic_table_wdg->setModel(model);
    m_atom_query_wdg->setModel(model);
}

void SetAtomWidget::updateCheckedButton()
{
    auto model = getModel();
    auto draw_tool = model->getDrawTool();
    auto atom_tool = model->getAtomTool();
    auto element = model->getElement();
    QAbstractButton* button = nullptr;

    if (draw_tool == DrawTool::ATOM) {
        if (atom_tool == AtomTool::QUERY) {
            auto query_int = model->getValueInt(ModelKey::ATOM_QUERY);
            ui->atom_query_btn->setEnumItem(query_int);
            button = ui->atom_query_btn;
        } else {
            if (m_button_element_bimap.right.count(element) == 1) {
                button = m_button_element_bimap.right.at(element);
            } else {
                ui->last_picked_element_btn->setElement(element);
                button = ui->last_picked_element_btn;
            }
        }
    }

    check_button_or_uncheck_group(button, ui->atom_group);
}

void SetAtomWidget::onModelValuePinged(ModelKey key, QVariant value)
{
    if (key == ModelKey::ELEMENT) {
        auto element = value.value<Element>();
        if (m_button_element_bimap.right.count(element) == 0) {
            ui->last_picked_element_btn->setElement(element);
        }
    }
}

void SetAtomWidget::updateWidgetsEnabled()
{
    AbstractDrawToolWidget::updateWidgetsEnabled();
    auto model = getModel();
    bool enable = (!model->hasActiveSelection() || model->hasAtomSelection());
    setEnabled(enable);
}

std::unordered_set<QAbstractButton*> SetAtomWidget::getCheckableButtons()
{
    auto buttons = AbstractDrawToolWidget::getCheckableButtons();
    for (auto button : ui->atom_group->buttons()) {
        buttons.insert(button);
    }
    return buttons;
}

void SetAtomWidget::onAtomButtonClicked(QAbstractButton* button)
{
    auto model = getModel();
    std::unordered_map<ModelKey, QVariant> kv_pairs;
    AtomTool atom_tool;
    if (button == ui->atom_query_btn) {
        atom_tool = AtomTool::QUERY;
        auto atom_query = AtomQuery(ui->atom_query_btn->getEnumItem());
        kv_pairs.emplace(ModelKey::ATOM_QUERY, QVariant::fromValue(atom_query));
    } else {
        atom_tool = AtomTool::ELEMENT;
        Element element = getElementForButton(button);
        kv_pairs.emplace(ModelKey::ELEMENT, QVariant::fromValue(element));
    }

    if (model->hasActiveSelection()) {
        // Do not alter the model, but ping the desired tool to indicate that we
        // want to replace the selection
        auto pair = *kv_pairs.begin();
        model->pingValue(pair.first, pair.second);
    } else {
        // No selection, so fully update the model
        kv_pairs.emplace(ModelKey::DRAW_TOOL,
                         QVariant::fromValue(DrawTool::ATOM));
        kv_pairs.emplace(ModelKey::ATOM_TOOL, QVariant::fromValue(atom_tool));

        model->setValues(kv_pairs);
    }
}

Element SetAtomWidget::getElementForButton(QAbstractButton* button) const
{
    if (button == ui->last_picked_element_btn) {
        return ui->last_picked_element_btn->getElement();
    } else {
        return m_button_element_bimap.left.at(button);
    }
}

SetAtomMenuWidget::SetAtomMenuWidget(QWidget* parent) : SetAtomWidget(parent)
{
    // Expand the margins so the buttons don't overlap the menu border
    ui->gridLayout->setContentsMargins(1, 1, 1, 1);
    // Hide the query types, given they are irrelevant for the menu widget
    ui->atom_query_btn->setVisible(false);

    // Pass signal up so that the parent menu can close on click
    connect(ui->atom_group,
            static_cast<void (QButtonGroup::*)(QAbstractButton*)>(
                &QButtonGroup::buttonClicked),
            this, &SetAtomMenuWidget::anyButtonClicked);
    connect(m_periodic_table_wdg, &PeriodicTableWidget::elementSelected, this,
            &SetAtomMenuWidget::anyButtonClicked);
}

SetAtomMenuWidget::~SetAtomMenuWidget() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/set_atom_widget.moc"
