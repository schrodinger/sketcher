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
 * Clear the originating button's hover and pressed state when an asynchronous
 * popup hides, preserving its checked state. This handles dismissal paths that
 * do not deliver a mouse leave event to the button. The popup owns this filter;
 * the button is tracked without taking ownership and may be destroyed first.
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
        if (event->type() == QEvent::Hide && m_button) {
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

void clear_button_state_on_popup_hide(QAbstractButton* button, QWidget* popup)
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
