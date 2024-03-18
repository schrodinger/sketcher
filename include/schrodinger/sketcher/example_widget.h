// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#include <memory>

#include <QWidget>

#include "schrodinger/sketcher/definitions.h"

namespace Ui
{
class ExampleWidgetForm;
}

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ExampleWidget : public QWidget
{
    Q_OBJECT
  public:
    ExampleWidget(QWidget* parent = nullptr);
    ~ExampleWidget();

  protected:
    std::unique_ptr<Ui::ExampleWidgetForm> m_ui;
};

} // namespace sketcher
} // namespace schrodinger
