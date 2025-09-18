#include "schrodinger/sketcher/widget/select_options_widget.h"

#include <functional>

#include "schrodinger/sketcher/widget/widget_utils.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_select_options_widget.h"
#include "schrodinger/sketcher/widget/selection_tool_popup.h"

namespace schrodinger
{
namespace sketcher
{

SelectOptionsWidget::SelectOptionsWidget(QWidget* parent) : SketcherView(parent)
{
    ui.reset(new Ui::SelectOptionsWidget());
    ui->setupUi(this);

    ui->title_lbl->setStyleSheet(PALETTE_TITLE_STYLE);
    for (auto btn : {ui->select_all_btn, ui->clear_selection_btn,
                     ui->invert_selection_btn}) {
        btn->setStyleSheet(TEXT_LINK_STYLE);
    }
    m_selection_tool1_widget = new SelectionToolPopup(this);
    m_selection_tool2_widget = new SelectionToolPopup(this);

    ui->select_tool_btn_1->setPopupWidget(m_selection_tool1_widget);
    ui->select_tool_btn_2->setPopupWidget(m_selection_tool2_widget);

    ui->select_tool_btn_1->setEnumItem(
        static_cast<int>(SelectionTool::RECTANGLE));
    ui->select_tool_btn_2->setEnumItem(static_cast<int>(SelectionTool::LASSO));

    connect(ui->select_tool_group,
            static_cast<void (QButtonGroup::*)(QAbstractButton*)>(
                &QButtonGroup::buttonClicked),
            this, &SelectOptionsWidget::onSelectButtonClicked);
    connect(ui->select_all_btn, &QToolButton::clicked, this,
            &SelectOptionsWidget::selectAllRequested);
    connect(ui->clear_selection_btn, &QToolButton::clicked, this,
            &SelectOptionsWidget::clearSelectionRequested);
    connect(ui->invert_selection_btn, &QToolButton::clicked, this,
            &SelectOptionsWidget::invertSelectionRequested);
    connect(ui->move_rotate_btn, &QToolButton::clicked, this,
            &SelectOptionsWidget::onMoveButtonClicked);
    connect(ui->erase_btn, &QToolButton::clicked, this,
            &SelectOptionsWidget::onEraseButtonClicked);
}

void SelectOptionsWidget::setModel(SketcherModel* model)
{
    SketcherView::setModel(model);
    m_selection_tool1_widget->setModel(model);
    m_selection_tool2_widget->setModel(model);
}

SelectOptionsWidget::~SelectOptionsWidget()
{
    delete m_selection_tool1_widget;
    delete m_selection_tool2_widget;
}

void SelectOptionsWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    auto has_contents = !model->sceneIsEmpty();
    auto has_selection = model->hasActiveSelection();

    // Something must be drawn for these buttons to do anything
    ui->select_tool_btn_1->setEnabled(has_contents);
    ui->select_tool_btn_2->setEnabled(has_contents);
    ui->erase_btn->setEnabled(has_contents);

    // Move tool should only be enabled if there are atoms or non-molecular
    // objects selected
    bool movable_selection = !has_selection || model->hasAtomSelection() ||
                             model->hasNonMolecularObjectSelection();
    ui->move_rotate_btn->setEnabled(has_contents && movable_selection);

    ui->select_all_btn->setEnabled(has_contents && !model->allItemsSelected());
    ui->clear_selection_btn->setEnabled(has_selection);
    ui->invert_selection_btn->setEnabled(has_selection);

    // Update background color of move/erase tools
    ui->move_frame->setStyleSheet(has_selection ? SELECTION_ACTIVE_STYLE : "");
}

void SelectOptionsWidget::updateCheckState()
{
    auto model = getModel();
    auto draw_tool = model->getDrawTool();

    QAbstractButton* select_button = nullptr;
    QAbstractButton* move_button = nullptr;
    QAbstractButton* erase_button = nullptr;
    if (draw_tool == DrawTool::SELECT) {
        auto button_id = model->getValueInt(ModelKey::SELECTION_TOOL);
        // check each modular button in the group to find if one is currently
        // set to this ID
        for (auto* btn : ui->select_tool_group->buttons()) {
            auto mod_button = dynamic_cast<ModularToolButton*>(btn);
            if (mod_button != nullptr &&
                mod_button->getEnumItem() == button_id) {
                select_button = btn;
                break;
            }
        }
    } else if (draw_tool == DrawTool::MOVE_ROTATE) {
        move_button = ui->move_rotate_btn;
    } else if (draw_tool == DrawTool::ERASE) {
        erase_button = ui->erase_btn;
    }
    check_button_or_uncheck_group(select_button, ui->select_tool_group);
    check_button_or_uncheck_group(move_button, ui->move_tool_group);
    check_button_or_uncheck_group(erase_button, ui->erase_tool_group);
}

void SelectOptionsWidget::onSelectButtonClicked(QAbstractButton* button)
{
    auto button_cast = dynamic_cast<ModularToolButton*>(button);
    if (!button_cast) {
        throw std::runtime_error("Invalid button type");
    }
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::SELECT)},
        {ModelKey::SELECTION_TOOL,
         QVariant::fromValue(SelectionTool(button_cast->getEnumItem()))},
    };
    getModel()->setValues(kv_pairs);
}

void SelectOptionsWidget::onMoveButtonClicked()
{
    getModel()->setValue(ModelKey::DRAW_TOOL, DrawTool::MOVE_ROTATE);
}

void SelectOptionsWidget::onEraseButtonClicked()
{
    auto model = getModel();
    if (model->hasActiveSelection()) {
        // do not change the model
        model->pingValue(ModelKey::DRAW_TOOL, DrawTool::ERASE);
    } else {
        // update the model
        model->setValue(ModelKey::DRAW_TOOL, DrawTool::ERASE);
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/select_options_widget.moc"
