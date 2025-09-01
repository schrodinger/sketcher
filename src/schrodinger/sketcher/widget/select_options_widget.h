#pragma once

#include <memory>

#include <QWidget>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"
#include <QAbstractButton>

namespace Ui
{
class SelectOptionsWidget;
}

namespace schrodinger
{
namespace sketcher
{

class SelectionToolPopup;

class SKETCHER_API SelectOptionsWidget : public SketcherView
{
    Q_OBJECT

  public:
    SelectOptionsWidget(QWidget* parent = nullptr);
    ~SelectOptionsWidget();

    void setModel(SketcherModel* model) override;

  signals:
    void selectAllRequested();
    void clearSelectionRequested();
    void invertSelectionRequested();

  protected:
    std::unique_ptr<Ui::SelectOptionsWidget> ui;
    SelectionToolPopup* m_selection_tool1_widget = nullptr;
    SelectionToolPopup* m_selection_tool2_widget = nullptr;

  protected slots:
    void updateWidgetsEnabled() override;

    void onSelectButtonClicked(QAbstractButton* button);
    void onMoveButtonClicked();
    void onEraseButtonClicked();

  private:
    void updateCheckState() override;
};

} // namespace sketcher
} // namespace schrodinger
