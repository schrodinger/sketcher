// @copyright Schrodinger, LLC - All Rights Reserved

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

class SKETCHER_API SketcherWidget : public QWidget
{
    Q_OBJECT
  public:
    SketcherWidget(QWidget* parent = nullptr);
    ~SketcherWidget();

  protected:
    std::unique_ptr<Ui::SketcherWidgetForm> m_ui;
};

} // namespace sketcher
} // namespace schrodinger
