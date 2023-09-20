#include "schrodinger/sketcher/sketcher_widget.h"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <QGraphicsPixmapItem>
#include <QWidget>
#include <fmt/format.h>

#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);

    // Hello World using RDKit, maeparser, fmt
    auto mol = RDKit::SmilesToMol("CC");
    auto mae_block = RDKit::MaeWriter::getText(*mol);
    RDKit::MaeMolSupplier reader;
    reader.setData(mae_block, false, false);
    auto natoms = reader.next()->getNumAtoms();
    setWindowTitle(fmt::format("Hello World: {}", natoms).c_str());
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"
