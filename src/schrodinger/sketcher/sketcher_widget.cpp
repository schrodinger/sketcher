#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui_sketcher_widget.h"
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
    ui->select_options_wdg->setModel(m_sketcher_model);
    ui->draw_tools_wdg->setModel(m_sketcher_model);
    ui->fragment_wdg->setModel(m_sketcher_model);

    setStyleSheet("QToolButton:checked {"
                  "background-color: rgb(190, 190, 190);"
                  "}");
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger
