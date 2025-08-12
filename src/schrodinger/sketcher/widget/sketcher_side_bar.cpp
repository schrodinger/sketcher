#include "schrodinger/sketcher/widget/sketcher_side_bar.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_side_bar.h"

namespace schrodinger
{
namespace sketcher
{

SketcherSideBar::SketcherSideBar(QWidget* parent) : SketcherView(parent)
{
    ui.reset(new Ui::SketcherSideBar());
    ui->setupUi(this);

    // Force maximum width
    setMaximumWidth(100);

    // Forward required slots
    connect(ui->select_options_wdg, &SelectOptionsWidget::selectAllRequested,
            this, &SketcherSideBar::selectAllRequested);
    connect(ui->select_options_wdg,
            &SelectOptionsWidget::clearSelectionRequested, this,
            &SketcherSideBar::clearSelectionRequested);
    connect(ui->select_options_wdg,
            &SelectOptionsWidget::invertSelectionRequested, this,
            &SketcherSideBar::invertSelectionRequested);
}

SketcherSideBar::~SketcherSideBar() = default;

void SketcherSideBar::setModel(SketcherModel* model)
{
    ui->select_options_wdg->setModel(model);
    ui->draw_tools_wdg->setModel(model);
    ui->ring_tool_wdg->setModel(model);
    ui->enumeration_tool_wdg->setModel(model);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/sketcher_side_bar.moc"
