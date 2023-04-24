#include "schrodinger/sketcher/widget/select_options_widget.h"

#include <functional>

#include "schrodinger/sketcher/qt_utils.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_select_options_widget.h"

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

    // Init buttons
    ui->select_tool_group->setId(ui->select_square_btn,
                                 static_cast<int>(SelectionTool::RECTANGLE));
    ui->select_tool_group->setId(ui->select_lasso_btn,
                                 static_cast<int>(SelectionTool::LASSO));
    connect(ui->select_tool_group, &QButtonGroup::idClicked, this,
            &SelectOptionsWidget::onSelectButtonClicked);
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

SelectOptionsWidget::~SelectOptionsWidget() = default;

void SelectOptionsWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    auto has_contents = !model->sceneIsEmpty();
    auto has_selection = model->hasActiveSelection();

    // Something must be drawn for these buttons to do anything
    ui->select_square_btn->setEnabled(has_contents);
    ui->select_lasso_btn->setEnabled(has_contents);
    ui->erase_btn->setEnabled(has_contents);

    // Move tool should only enabled if there are atoms selected
    bool movable_selection = !has_selection || model->hasAtomSelection();
    ui->move_rotate_btn->setEnabled(has_contents && movable_selection);

    ui->select_all_btn->setEnabled(has_contents);
    ui->clear_selection_btn->setEnabled(has_selection);
    ui->invert_selection_btn->setEnabled(has_selection);

    // Update background color of move/erase tools
    ui->move_frame->setStyleSheet(has_selection ? SELECTION_ACTIVE_STYLE : "");
}

void SelectOptionsWidget::updateCheckState()
{
    auto model = getModel();
    auto draw_tool = model->getDrawTool();

    QAbstractButton* button = nullptr;
    if (draw_tool == DrawTool::SELECT) {
        auto button_id = model->getValueInt(ModelKey::SELECTION_TOOL);
        button = ui->select_tool_group->button(button_id);
    }
    check_button_or_uncheck_group(button, ui->select_tool_group);
    ui->move_rotate_btn->setChecked(draw_tool == DrawTool::MOVE_ROTATE);
    ui->erase_btn->setChecked(draw_tool == DrawTool::ERASE);
}

void SelectOptionsWidget::onSelectButtonClicked(int button_id)
{
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::SELECT)},
        {ModelKey::SELECTION_TOOL,
         QVariant::fromValue(SelectionTool(button_id))},
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
