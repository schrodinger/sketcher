#pragma once

#include <memory>
#include <QWidget>
#include "schrodinger/sketcher/definitions.h"

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
    ~SketcherWidget();

  private:
    std::unique_ptr<Ui::SketcherWidgetForm> ui;
    SketcherModel* m_sketcher_model = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
