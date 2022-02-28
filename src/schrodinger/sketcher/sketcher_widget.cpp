#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"
#include "schrodinger/sketcher/sketcher_css_style.h" // TOOL_BUTTON_STYLE
#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    setProperty("SKIP_LEGACY_STYLING", true);
    ui.reset(new Ui::SketcherWidgetForm());
    ui->setupUi(this);

    m_sketcher_model = new SketcherModel(this);
    ui->top_bar_wdg->setModel(m_sketcher_model);
    ui->side_bar_wdg->setModel(m_sketcher_model);

    setStyleSheet(TOOL_BUTTON_STYLE);
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger
