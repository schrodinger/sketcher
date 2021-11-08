#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui_sketcher_widget.h"
#include "schrodinger/sketcher/sketcher_model.h"

using namespace schrodinger::sketcher;

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    setProperty("SKIP_LEGACY_STYLING", true);
    ui = new Ui::SketcherWidgetForm();
    ui->setupUi(this);

    m_sketcher_model = new SketcherModel(this);
    ui->top_bar_wdg->setModel(m_sketcher_model);
}

SketcherWidget::~SketcherWidget()
{
    delete ui;
}
