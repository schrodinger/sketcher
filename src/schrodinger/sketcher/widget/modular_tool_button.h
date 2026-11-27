#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"

namespace schrodinger
{
namespace sketcher
{

class ModularPopup;

/**
 * A button that allows you to swap out its appearance/functionality via popup.
 *
 * Meant to be used with the `ModularPopup` class. Additionally, when a new
 * functionality is selected via the popup, this button will replace its icon
 * with that of the selected button, as well as its assigned enum-based integer
 * value.
 */
class SKETCHER_API ModularToolButton : public ToolButtonWithPopup
{
    Q_OBJECT
  public:
    ModularToolButton(QWidget* parent = nullptr);
    void setPopupWidget(QWidget* popup) override;

    /**
     * @return The integer enum value currently assigned to this button.
     */
    int getEnumItem();

    /**
     * Assign the integer enum value for this button.
     *
     * The integer enum value on this button can be used to determine its
     * intended functionality.
     */
    void setEnumItem(int enum_int);

  signals:
    void enumChanged();

  protected:
    /**
     * Update button text and icon to match corresponding button in popup.
     */
    void updateButton();
    int m_enum_int = -1;

  protected slots:
    /**
     * Respond to a change in the popup widget's selection state by updating the
     * enum integer and icon assigned to this button.
     */
    void onPopupSelectionChanged(int enum_int);
};

} // namespace sketcher
} // namespace schrodinger
