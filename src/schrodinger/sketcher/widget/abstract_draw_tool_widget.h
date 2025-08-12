#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

class QAbstractButton;

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API AbstractDrawToolWidget : public SketcherView
{
    Q_OBJECT

  public:
    AbstractDrawToolWidget(QWidget* parent = nullptr);
    ~AbstractDrawToolWidget() = default;

  protected:
    /**
     * Update whether draw tool buttons are checkable.
     */
    void updateCheckState() override;

    /**
     * Assign the appropriate check state based on the state of the model.
     */
    virtual void updateCheckedButton() = 0;

    /**
     * @return All buttons that should be checkable in draw mode.
     */
    virtual std::unordered_set<QAbstractButton*> getCheckableButtons();
};

} // namespace sketcher
} // namespace schrodinger
