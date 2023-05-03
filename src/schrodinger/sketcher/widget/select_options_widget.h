#pragma once

#include <memory>

#include <QWidget>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

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

  signals:
    void selectAllRequested();
    void clearSelectionRequested();
    void invertSelectionRequested();

  protected:
    std::unique_ptr<Ui::SelectOptionsWidget> ui;

  protected slots:
    void updateWidgetsEnabled() override;

    void onSelectButtonClicked(int button_id);
    void onMoveButtonClicked();
    void onEraseButtonClicked();

  private:
    void updateCheckState() override;
};

} // namespace sketcher
} // namespace schrodinger
