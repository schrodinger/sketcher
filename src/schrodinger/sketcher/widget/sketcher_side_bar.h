#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

class QAbstractButton;
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
    void updateWidgetsEnabled() override;
    void updateCheckState() override;

    /**
     * Disconnect signals from each widget's updateWidgetsEnabled()
     * slot; called when selection-only mode is enabled in SketcherWidget.
     */
    void disconnectAllUpdateWidgetsEnabled();

  signals:
    void selectAllRequested();
    void clearSelectionRequested();
    void invertSelectionRequested();

  protected:
    std::unique_ptr<Ui::SketcherSideBar> ui;

    void onAtomisticOrMonomerButtonClicked(QAbstractButton* button);

    DrawTool m_previous_atomistic_draw_tool = DrawTool::ATOM;
};

} // namespace sketcher
} // namespace schrodinger
