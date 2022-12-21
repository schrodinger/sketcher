#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

class QWidget;

namespace Ui
{
class SketcherSideBar;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Side bar containing palettes for the sketcher widget
 */
class SKETCHER_API SketcherSideBar : public SketcherView
{
    Q_OBJECT

  public:
    SketcherSideBar(QWidget* parent = nullptr);
    ~SketcherSideBar();

    void setModel(SketcherModel* model) override;

  signals:
    void selectAllRequested();
    void clearSelectionRequested();
    void invertSelectionRequested();

  protected:
    std::unique_ptr<Ui::SketcherSideBar> ui;
};

} // namespace sketcher
} // namespace schrodinger
