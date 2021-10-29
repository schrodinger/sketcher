#pragma once
#include "schrodinger/sketcher/definitions.h"
#include <QWidget>

namespace Ui
{
class SketcherWidgetForm;
}

namespace schrodinger
{
namespace sketcher
{

class SketcherModel;

/**
 * Sketcher widget meant for use in LiveDesign and Maestro.
 */
class SKETCHER_API SketcherWidget : public QWidget
{
  public:
    SketcherWidget(QWidget* parent = nullptr);

  private:
    Ui::SketcherWidgetForm* ui = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
