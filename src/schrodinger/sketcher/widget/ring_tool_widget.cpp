#include "schrodinger/sketcher/widget/ring_tool_widget.h"

#include <functional>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_ring_tool_widget.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

RingToolWidget::RingToolWidget(QWidget* parent) : SketcherView(parent)
{
    ui.reset(new Ui::RingToolWidget());
    ui->setupUi(this);

    ui->group->setId(ui->cyclopropane_btn,
                     static_cast<int>(RingTool::CYCLOPROPANE));
    ui->group->setId(ui->cyclobutane_btn,
                     static_cast<int>(RingTool::CYCLOBUTANE));
    ui->group->setId(ui->cyclopentane_btn,
                     static_cast<int>(RingTool::CYCLOPENTANE));
    ui->group->setId(ui->cyclopentadiene_btn,
                     static_cast<int>(RingTool::CYCLOPENTADIENE));
    ui->group->setId(ui->cyclohexane_btn,
                     static_cast<int>(RingTool::CYCLOHEXANE));
    ui->group->setId(ui->benzene_btn, static_cast<int>(RingTool::BENZENE));
    ui->group->setId(ui->cycloheptane_btn,
                     static_cast<int>(RingTool::CYCLOHEPTANE));
    ui->group->setId(ui->cyclooctane_btn,
                     static_cast<int>(RingTool::CYCLOOCTANE));
    connect(ui->group, &QButtonGroup::idClicked, this,
            &RingToolWidget::onButtonClicked);
}

RingToolWidget::~RingToolWidget() = default;

void RingToolWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    setEnabled(!model->hasActiveSelection());
}

void RingToolWidget::updateCheckState()
{
    auto model = getModel();
    auto draw_tool = model->getDrawTool();
    QAbstractButton* button = nullptr;
    if (draw_tool == DrawTool::RING) {
        auto button_id = model->getValueInt(ModelKey::RING_TOOL);
        button = ui->group->button(button_id);
    }

    check_button_or_uncheck_group(button, ui->group);
}

void RingToolWidget::onButtonClicked(int button_id)
{
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::RING)},
        {ModelKey::RING_TOOL, QVariant::fromValue(RingTool(button_id))},
    };
    getModel()->setValues(kv_pairs);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/ring_tool_widget.moc"
