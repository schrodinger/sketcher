#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"

#include <QPainter>
#include <QPaintEvent>
#include <QStyleOption>

namespace schrodinger
{
namespace sketcher
{

SketcherViewWithWasmOutlineFix::SketcherViewWithWasmOutlineFix(
    QWidget* parent) :
    SketcherView(parent)
{
}

SketcherViewWithWasmOutlineFix::~SketcherViewWithWasmOutlineFix() = default;

void SketcherViewWithWasmOutlineFix::paintEvent(QPaintEvent*)
{
    QStyleOption opt;
    opt.initFrom(this);
    QPainter p(this);
    style()->drawPrimitive(QStyle::PE_Widget, &opt, &p, this);
}

} // namespace sketcher
} // namespace schrodinger
