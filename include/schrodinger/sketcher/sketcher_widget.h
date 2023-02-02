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

class Scene;
class SketcherModel;

/**
 * Sketcher widget meant for use in LiveDesign and Maestro.
 */
class SKETCHER_API SketcherWidget : public QWidget
{
    Q_OBJECT

  public:
    SketcherWidget(QWidget* parent = nullptr);
    ~SketcherWidget();

  private:
    std::unique_ptr<Ui::SketcherWidgetForm> m_ui;
    Scene* m_scene = nullptr;
    SketcherModel* m_sketcher_model = nullptr;

    /**
     *  Connects slots to the model and various widget tool bars
     */
    void connectTopBarSlots();
    void connectSideBarSlots();
};

} // namespace sketcher
} // namespace schrodinger
