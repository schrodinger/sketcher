#include "schrodinger/sketcher/widget/widget_utils.h"

#include <QAbstractButton>
#include <QApplication>
#include <QButtonGroup>
#include <QEvent>
#include <QMouseEvent>
#include <QPointer>

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Manage dismissal of a button's asynchronous popup: suppress replay of presses
 * over the button and clear hover and pressed state when the popup hides. The
 * popup owns this filter; the button is tracked without taking ownership and
 * may be destroyed first.
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
        if (event->type() == QEvent::MouseButtonPress && m_button) {
            auto* mouse_event = static_cast<QMouseEvent*>(event);
            if (QRect(m_button->mapToGlobal(QPoint()), m_button->size())
                    .contains(mouse_event->globalPosition().toPoint())) {
                // Match Qt's built-in button menus: do not replay a dismissing
                // press onto the button, since that will reopen the popup as
                // soon as the button is released. Note that Qt will clear the
                // WA_NoMouseReplay attribute when the popup closes, so we don't
                // have to worry about clearing it ourselves.
                static_cast<QWidget*>(watched)->setAttribute(
                    Qt::WA_NoMouseReplay);
            }
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
