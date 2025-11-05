#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"

class QPushButton;

namespace Ui
{
class PeriodicTableForm;
}

namespace schrodinger
{
namespace sketcher
{

enum class Element;
class SketcherModel;

/**
 * A free-floating widget that contains a button corresponding to each element
 * on the periodic table.
 */
class SKETCHER_API PeriodicTableWidget : public SketcherViewWithWasmOutlineFix
{
    Q_OBJECT
  public:
    PeriodicTableWidget(QWidget* parent = nullptr);
    ~PeriodicTableWidget();

    /**
     * @param close_on_click whether this widget should close when the user
     * clicks any button.
     */
    void setCloseOnClick(bool close_on_click);

  signals:
    /**
     * Emitted when the user clicks an element button.
     */
    void elementSelected(Element element);

  protected slots:
    void onButtonClicked(int button_id);

  protected:
    std::unique_ptr<Ui::PeriodicTableForm> ui;

  private:
    /**
     * Whether to close the widget when the user clicks any element button.
     */
    bool m_close_on_click = true;
};

} // namespace sketcher
} // namespace schrodinger
