#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

namespace Ui
{
class RingToolWidget;
}

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API RingToolWidget : public SketcherView
{
    Q_OBJECT

  public:
    RingToolWidget(QWidget* parent = nullptr);
    ~RingToolWidget();

  protected:
    std::unique_ptr<Ui::RingToolWidget> ui;
    void updateWidgetsEnabled() override;

  protected slots:
    /**
     * Respond to the user interacting with a ring button by updating the model.
     *
     * @param button_id The ID of the desired ring button; this should
     * correspond to the integer-cast value of the corresponding ring tool enum
     */
    void onButtonClicked(int button_id);

  private:
    void updateCheckState() override;
};

} // namespace sketcher
} // namespace schrodinger
