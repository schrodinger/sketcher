// @copyright Schrodinger, LLC - All Rights Reserved

#include "schrodinger/sketcher/sketcher_widget.h"

#include <QWidget>

#include "schrodinger/rdkit_extensions/example.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);
    if (rdkit_extensions::dependency_test("c1ccccc1")) {
        m_ui->label->setText(QString::fromStdString("Dependencies work!"));
    }
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"
