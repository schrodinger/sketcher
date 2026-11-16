#include "schrodinger/sketcher/widget/widget_utils.h"

#include <QAbstractButton>
#include <QApplication>
#include <QButtonGroup>
#include <QEvent>
#include <QPointer>

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Manage dismissal of a button's asynchronous popup: suppress replay of
 * dismissing presses and clear hover and pressed state when the popup hides.
 * The popup owns this filter; the button is tracked without taking ownership
 * and may be destroyed first.
 */
class PopupButtonStateFilter : public QObject
{
  public:
    PopupButtonStateFilter(QAbstractButton* button, QWidget* popup) :
        QObject(popup),
        m_button(button)
    {
        popup->installEventFilter(this);
    }

  protected:
    bool eventFilter(QObject* watched, QEvent* event) override
    {
        if (event->type() == QEvent::MouseButtonPress) {
            // A press that dismisses the popup should only dismiss it. By
            // default, Qt replays the press onto the widget underneath on
            // Windows and Linux, but not on macOS (where the platform
            // disables replay) or WASM (where popups can't grab the
            // keyboard, which also disables replay). Without this, a click on
            // the button would reopen the popup on release, and a click on the
            // canvas would also, e.g., add a monomer. Qt clears this attribute
            // before delivering each mouse event, so it only affects this
            // press, and Qt only replays presses outside the popup, so this
            // has no effect on presses inside of it.
            static_cast<QWidget*>(watched)->setAttribute(Qt::WA_NoMouseReplay);
        } else if (event->type() == QEvent::Hide && m_button) {
            // A popup grabs the mouse, so the originating button can miss its
            // leave event when the user dismisses it outside the client area.
            QEvent leave(QEvent::Leave);
            QApplication::sendEvent(m_button, &leave);
            m_button->setDown(false);
            m_button->update();
        }
        return QObject::eventFilter(watched, event);
    }

  private:
    QPointer<QAbstractButton> m_button;
};

} // namespace

void configure_button_popup(QAbstractButton* button, QWidget* popup)
{
    new PopupButtonStateFilter(button, popup);
}

/**
 * Uncheck all buttons in a button group.
 *
 * @param group The button group to uncheckify
 */
static void uncheck_all(QButtonGroup* group)
{
    bool exclusive = group->exclusive();
    group->setExclusive(false);
    for (auto button : group->buttons()) {
        button->setChecked(false);
    }
    group->setExclusive(exclusive);
}

void check_button_or_uncheck_group(QAbstractButton* button, QButtonGroup* group)
{
    if (button == nullptr) {
        uncheck_all(group);
    } else {
        button->setChecked(true);
    }
}

} // namespace sketcher
} // namespace schrodinger
