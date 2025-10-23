#include "schrodinger/sketcher/widget/draw_tools_widget.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_draw_tools_widget.h"
#include "schrodinger/sketcher/widget/bond_order_popup.h"
#include "schrodinger/sketcher/widget/bond_query_popup.h"
#include "schrodinger/sketcher/widget/stereo_bond_popup.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

using schrodinger::sketcher::BondTool;

namespace schrodinger
{
namespace sketcher
{

DrawToolsWidget::DrawToolsWidget(QWidget* parent) :
    AbstractDrawToolWidget(parent)
{
    ui.reset(new Ui::DrawToolsWidget());
    ui->setupUi(this);

    m_bond_order_wdg = new BondOrderPopup(this);
    m_bond_query_wdg = new BondQueryPopup(this);
    m_stereo_bond1_wdg = new StereoBondPopup(this);
    m_stereo_bond2_wdg = new StereoBondPopup(this);

    ui->bond_order_btn->setEnumItem(static_cast<int>(BondTool::DOUBLE));
    ui->bond_query_btn->setEnumItem(static_cast<int>(BondTool::AROMATIC));
    ui->stereo_bond1_btn->setEnumItem(static_cast<int>(BondTool::SINGLE_UP));
    ui->stereo_bond2_btn->setEnumItem(static_cast<int>(BondTool::SINGLE_DOWN));

    ui->bond_order_btn->setPopupWidget(m_bond_order_wdg);
    ui->bond_query_btn->setPopupWidget(m_bond_query_wdg);
    ui->stereo_bond1_btn->setPopupWidget(m_stereo_bond1_wdg);
    ui->stereo_bond2_btn->setPopupWidget(m_stereo_bond2_wdg);

    ui->bond_query_btn->setStyleSheet(BOND_QUERY_STYLE);

    // Not all of these IDs will be used, but because we need to assign some of
    // them manually we assign all of them manually in order to ensure that each
    // ID is unique
    ui->bond_group->setId(ui->single_bond_btn,
                          static_cast<int>(BondTool::SINGLE));
    ui->bond_group->setId(ui->bond_order_btn,
                          static_cast<int>(BondTool::DOUBLE));
    ui->bond_group->setId(ui->bond_query_btn,
                          static_cast<int>(BondTool::AROMATIC));
    ui->bond_group->setId(ui->stereo_bond1_btn,
                          static_cast<int>(BondTool::SINGLE_UP));
    ui->bond_group->setId(ui->stereo_bond2_btn,
                          static_cast<int>(BondTool::SINGLE_DOWN));
    ui->bond_group->setId(ui->atom_chain_btn,
                          static_cast<int>(BondTool::ATOM_CHAIN));

    ui->charge_group->setId(ui->decrease_charge_btn,
                            static_cast<int>(ChargeTool::DECREASE));
    ui->charge_group->setId(ui->increase_charge_btn,
                            static_cast<int>(ChargeTool::INCREASE));
}

DrawToolsWidget::~DrawToolsWidget()
{
    delete m_bond_order_wdg;
    delete m_bond_query_wdg;
    delete m_stereo_bond1_wdg;
    delete m_stereo_bond2_wdg;
}

void DrawToolsWidget::setModel(SketcherModel* model)
{
    AbstractDrawToolWidget::setModel(model);
    ui->set_atom_wdg->setModel(model);
    m_bond_order_wdg->setModel(model);
    m_bond_query_wdg->setModel(model);
    m_stereo_bond1_wdg->setModel(model);
    m_stereo_bond2_wdg->setModel(model);
}

void DrawToolsWidget::connectLocalSlots()
{
    AbstractDrawToolWidget::connectLocalSlots();

    connect(ui->bond_group,
            static_cast<void (QButtonGroup::*)(QAbstractButton*)>(
                &QButtonGroup::buttonClicked),
            this, &DrawToolsWidget::onBondButtonClicked);
    connect(ui->charge_group, &QButtonGroup::idClicked, this,
            &DrawToolsWidget::onChargeButtonClicked);
    connect(ui->explicit_h_btn, &QToolButton::clicked, this,
            &DrawToolsWidget::onExplicitHButtonClicked);
}

void DrawToolsWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    auto has_selection = model->hasActiveSelection();

    // Atom tools
    auto sel_has_atom = model->hasAtomSelection();
    bool enable_atom = (!has_selection || sel_has_atom);
    std::vector<QWidget*> widgets = {
        ui->increase_charge_btn, ui->decrease_charge_btn, ui->explicit_h_btn};
    for (auto wdg : widgets) {
        wdg->setEnabled(enable_atom);
    }

    // Bond tools
    bool sel_has_bond = model->hasBondSelection();
    bool enable_bond = (!has_selection || sel_has_bond);
    widgets = {ui->single_bond_btn, ui->bond_order_btn, ui->bond_query_btn,
               ui->stereo_bond1_btn, ui->stereo_bond2_btn};
    for (auto wdg : widgets) {
        wdg->setEnabled(enable_bond);
    }

    // Atom chain never available as edit action
    ui->atom_chain_btn->setEnabled(!has_selection);

    // Update title label and background color based on selection
    ui->atom_frame->setStyleSheet(sel_has_atom ? SELECTION_ACTIVE_STYLE : "");
    ui->bond_frame->setStyleSheet(sel_has_bond ? SELECTION_ACTIVE_STYLE : "");
}

std::unordered_set<QAbstractButton*> DrawToolsWidget::getCheckableButtons()
{
    auto buttons = AbstractDrawToolWidget::getCheckableButtons();
    for (auto group : {ui->bond_group, ui->charge_group}) {
        for (auto button : group->buttons()) {
            buttons.insert(button);
        }
    }
    buttons.insert(ui->explicit_h_btn);
    return buttons;
}

void DrawToolsWidget::updateCheckedButton()
{
    auto model = getModel();
    if (model == nullptr) {
        return;
    }

    QAbstractButton* bond_button = nullptr;
    QAbstractButton* charge_button = nullptr;
    QAbstractButton* explicit_h_button = nullptr;
    auto draw_tool = model->getDrawTool();

    if (draw_tool == DrawTool::BOND) {
        auto bond_tool = model->getBondTool();
        bond_button = getBondButton(bond_tool);
    } else if (draw_tool == DrawTool::CHARGE) {
        auto charge_tool_int = model->getValueInt(ModelKey::CHARGE_TOOL);
        charge_button = ui->charge_group->button(charge_tool_int);
    } else if (draw_tool == DrawTool::EXPLICIT_H) {
        explicit_h_button = ui->explicit_h_btn;
    }

    check_button_or_uncheck_group(bond_button, ui->bond_group);
    check_button_or_uncheck_group(charge_button, ui->charge_group);
    check_button_or_uncheck_group(explicit_h_button, ui->explicit_h_group);
}

QAbstractButton* DrawToolsWidget::getBondButton(BondTool tool)
{
    // Try finding a button already associated with the specified tool
    auto group = ui->bond_group;
    auto tool_int = static_cast<int>(tool);
    for (auto button : group->buttons()) {
        auto mod_button = dynamic_cast<ModularToolButton*>(button);
        if ((mod_button == nullptr && group->id(button) == tool_int) ||
            (mod_button != nullptr && mod_button->getEnumItem() == tool_int)) {
            return button;
        }
    }

    // Try assigning the enum to an existing modular button (note that the
    // second stereo bond popup is excluded, as it is redundant)
    std::unordered_set<ModularPopup*> popups = {
        m_bond_order_wdg, m_bond_query_wdg, m_stereo_bond1_wdg};
    for (auto popup : popups) {
        for (auto button : group->buttons()) {
            if (popup->getButtonIDs().count(tool_int) == 1) {
                auto mod_button = dynamic_cast<ModularToolButton*>(button);
                if (mod_button != nullptr &&
                    mod_button->getPopupWidget() == popup) {
                    mod_button->setEnumItem(tool_int);
                    return button;
                }
            }
        }
    }

    // The above methods should not fail, but just in case
    return nullptr;
}

void DrawToolsWidget::onBondButtonClicked(QAbstractButton* button)
{
    auto model = getModel();
    auto button_cast = dynamic_cast<ModularToolButton*>(button);
    BondTool bond_tool;
    if (button_cast == nullptr) {
        bond_tool = BondTool(ui->bond_group->id(button));
    } else {
        bond_tool = BondTool(button_cast->getEnumItem());
    }

    if (model->hasActiveSelection()) {
        // do not change the model
        model->pingValue(ModelKey::BOND_TOOL, bond_tool);
    } else {
        // update the model
        std::unordered_map<ModelKey, QVariant> kv_pairs = {
            {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::BOND)},
            {ModelKey::BOND_TOOL, QVariant::fromValue(bond_tool)},
        };
        model->setValues(kv_pairs);
    }
}

void DrawToolsWidget::onChargeButtonClicked(int button_id)
{
    auto model = getModel();
    auto charge_tool = ChargeTool(button_id);
    if (model->hasActiveSelection()) {
        // do not change the model
        model->pingValue(ModelKey::CHARGE_TOOL, charge_tool);
    } else {
        // update the model
        std::unordered_map<ModelKey, QVariant> kv_pairs = {
            {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::CHARGE)},
            {ModelKey::CHARGE_TOOL, QVariant::fromValue(charge_tool)},
        };
        model->setValues(kv_pairs);
    }
}

void DrawToolsWidget::onExplicitHButtonClicked()
{
    auto model = getModel();
    if (model->hasActiveSelection()) {
        // do not change the model
        model->pingValue(ModelKey::DRAW_TOOL,
                         QVariant::fromValue(DrawTool::EXPLICIT_H));
    } else {
        // update the model
        model->setValue(ModelKey::DRAW_TOOL, DrawTool::EXPLICIT_H);
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/draw_tools_widget.moc"
