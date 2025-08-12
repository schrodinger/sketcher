#pragma once
#include <QToolButton>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ToolButtonWithPopup : public QToolButton
{
    Q_OBJECT
  public:
    ToolButtonWithPopup(QWidget* parent = nullptr);

    /**
     * Assign the popup widget for this button.
     *
     * This widget can only be set a single time.
     *
     * @throws std::runtime_error If this method is called when the popup widget
     * has already been set
     * @param popup_wdg The new popup widget to use for this button
     */
    virtual void setPopupWidget(QWidget* popup_wdg);

    /**
     * @return the popup widget assigned to this button
     */
    QWidget* getPopupWidget() const;

    /**
     * Set the duration the user must hold left click before the popup appears.
     *
     * @param popup_delay The new popup delay duration, in ms
     */
    void setPopupDelay(float popup_delay);

    /**
     * Show or hide the popup indicator (the wedge at the bottom right).
     *
     * @param show Whether to show or hide the indicator
     */
    void showPopupIndicator(bool show);

    /**
     * Apply additional style sheet parameters to the default style used for
     * this widget.
     *
     * @param text Text of the desired additional style sheet
     */
    void setStyleSheet(QString text);

  protected:
    QWidget* m_popup_wdg = nullptr;
    QString m_custom_style_sheet;

    /**
     * For tracking button check state.
     *
     * Used to determine whether the button has been clicked once it's already
     * been checked. This is impossible to do with built-in state, because once
     * the `clicked()` signal has been emitted, that check state has already
     * been updated.
     */
    bool m_was_checked_on_press = false;

    /**
     * Show the popup widget.
     *
     * By default, place the top left corner of the widget at the bottom left
     * corner of this button. If the widget would be rendered partially off-
     * screen, move it up and left as necessary. If it needs to be moved up,
     * show the popup fully above this button.
     */
    void showPopup();
    void startPopupTimer();

    /**
     * Apply the style sheet saved in `m_style_sheet`.
     *
     * If this is done during construction, the style sheet will be overwritten
     * on any instances that appear in a .ui file. Consequently, this method is
     * called from `setPopupWidget()`.
     */
    void updateStyle();

    /**
     * Override paintEvent() to ensure that a menu indicator is always drawn.
     *
     * Normally, the menu indicator would only be rendered if a menu is set on
     * the button. However, we don't use a conventional menu as our popup, so we
     * must trick the painter into thinking we have one.
     */
    void paintEvent(QPaintEvent* event) override;

    /**
     * Respond when the button is clicked by showing the popup if already
     * checked.
     *
     * @param checked Whether the button is checked
     */
    void onClicked(bool checked);

    /**
     * Update the state of `m_was_checked_on_press` when the user presses the
     * button.
     *
     * This must be done at this time (prior to a click) to correctly capture
     * the checkstate of the button before the click occurs. After the click
     * occurs, the button will always be checked.
     */
    void onPressed();

  private:
    QTimer* m_popup_timer = nullptr;
    float m_popup_delay = 250;
    bool m_show_corner_arrow = true;
};

} // namespace sketcher
} // namespace schrodinger
