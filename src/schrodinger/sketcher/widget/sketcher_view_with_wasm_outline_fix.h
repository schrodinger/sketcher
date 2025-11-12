#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A SketcherView that uses a custom paintEvent method to fix SKETCH-1553, which
 * causes popups to be drawn without an outline on WASM builds.
 */
class SKETCHER_API SketcherViewWithWasmOutlineFix : public SketcherView
{
  public:
    SketcherViewWithWasmOutlineFix(QWidget* parent = nullptr);
    ~SketcherViewWithWasmOutlineFix();
    void paintEvent(QPaintEvent*) override;
};

} // namespace sketcher
} // namespace schrodinger
