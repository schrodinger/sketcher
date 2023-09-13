#include "schrodinger/sketcher/sketcher_widget.h"

#include <GraphMol/FileParsers/MolSupplier.h>
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

    // Test using RDKit, maeparser, coordgen, fmt
    RDKit::MaeMolSupplier reader;
    auto text = R"(f_m_ct {
     s_m_title
     i_m_ct_format
     :::
     Structure1
      2
     m_atom[2] {
      # First column is atom index #
      i_m_mmod_type
      r_m_x_coord
      r_m_y_coord
      r_m_z_coord
      i_m_color
      i_m_atomic_number
      s_m_color_rgb
      i_rdk_index
      :::
      1 57 -0.076559 -0.161239 1.028628 9 17 006400  0
      2 57 0.076559 0.161239 -1.028628 9 17 006400  1
      :::
     }
     m_bond[1] {
      # First column is bond index #
      i_m_from
      i_m_to
      i_m_order
      :::
      1 1 2 1
      :::
     }
    })";
    reader.setData(text, false, false);
    auto mol = reader.next();
    // RDKit::CoordGen::addCoords(*mol);
    setWindowTitle(fmt::format("CC={} atoms", mol->getNumAtoms()).c_str());

    auto watermark = new QGraphicsPixmapItem();
    watermark->setPixmap(QPixmap(":icons/watermark.svg"));
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"