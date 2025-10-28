#include "schrodinger/sketcher/widget/sketcher_view_with_white_background.h"

#include <QPainter>
#include <QPaintEvent>
#include <QStyleOption>

namespace schrodinger
{
namespace sketcher
{

SketcherViewWithWhiteBackground::SketcherViewWithWhiteBackground(
    QWidget* parent) :
    SketcherView(parent)
{
}

SketcherViewWithWhiteBackground::~SketcherViewWithWhiteBackground() = default;

void SketcherViewWithWhiteBackground::paintEvent(QPaintEvent*)
{
    QStyleOption opt;
    opt.initFrom(this);
    QPainter p(this);
    style()->drawPrimitive(QStyle::PE_Widget, &opt, &p, this);
}

} // namespace sketcher
} // namespace schrodinger
