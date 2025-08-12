#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

namespace Ui
{
class EnumerationToolWidget;
}

namespace schrodinger
{
namespace sketcher
{

class ReactionPopup;
class RGroupPopup;

class SKETCHER_API EnumerationToolWidget : public SketcherView
{
    Q_OBJECT

  public:
    EnumerationToolWidget(QWidget* parent = nullptr);
    ~EnumerationToolWidget();
    void setModel(SketcherModel* model) override;

  protected:
    std::unique_ptr<Ui::EnumerationToolWidget> ui;
    ReactionPopup* m_reaction_popup = nullptr;
    RGroupPopup* m_rgroup_popup = nullptr;
    void updateWidgetsEnabled() override;

    /**
     * Update the active mode of modular buttons to match the model state.
     */
    void updateButtons();

  protected slots:
    void onRGroupButtonClicked();
    void onAttachmentPointButtonClicked();
    void onReactionButtonClicked();

  private:
    void updateCheckState() override;
};

} // namespace sketcher
} // namespace schrodinger
