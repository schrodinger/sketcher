#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A SketcherView that uses a custom paintEvent method to allow for a white
 * background when one is specified in the style sheet
 */
class SKETCHER_API SketcherViewWithWhiteBackground : public SketcherView
{
  public:
    SketcherViewWithWhiteBackground(QWidget* parent = nullptr);
    ~SketcherViewWithWhiteBackground();
    void paintEvent(QPaintEvent*) override;
};

} // namespace sketcher
} // namespace schrodinger
