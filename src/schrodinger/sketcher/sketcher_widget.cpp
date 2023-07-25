#include "schrodinger/sketcher/sketcher_widget.h"

// #include <GraphMol/SmilesParse/SmilesParse.h>
#include <QGraphicsPixmapItem>
#include <QWidget>
// #include <fmt/format.h>

#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);

    // Test using RDKit and fmt
    // auto mol = RDKit::SmilesToMol("c1ccccc1");
    // fmt::format("Benzene has {} atoms", mol.getNumAtoms());

    auto watermark = new QGraphicsPixmapItem();
    watermark->setPixmap(QPixmap(":icons/watermark.svg"));
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"