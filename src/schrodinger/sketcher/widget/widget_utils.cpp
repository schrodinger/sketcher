#include "schrodinger/sketcher/widget/widget_utils.h"

#include <algorithm>

#include <QAbstractButton>
#include <QApplication>
#include <QButtonGroup>
#include <QEvent>
#include <QMouseEvent>
#include <QPointer>
#include <QScreen>

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

QPoint clamp_position_to_bounds(const QRect& rect, const QRect& bounds)
{
    // Apply the bottom/right limit first so that the top/left limit wins when
    // the rectangle is too large to fit
    const int x =
        std::max(std::min(rect.left(), bounds.right() - rect.width() + 1),
                 bounds.left());
    const int y =
        std::max(std::min(rect.top(), bounds.bottom() - rect.height() + 1),
                 bounds.top());
    return {x, y};
}

QPoint calc_popup_position(const QSize& popup_size, const QRect& anchor_rect,
                           const QRect& bounds)
{
    const int room_below = bounds.bottom() - anchor_rect.bottom();
    const int room_above = anchor_rect.top() - bounds.top();
    const int below_y = anchor_rect.bottom() + 1;
    const int above_y = anchor_rect.top() - popup_size.height();
    int y;
    if (popup_size.height() <= room_below) {
        y = below_y;
    } else if (popup_size.height() <= room_above) {
        y = above_y;
    } else {
        y = room_below >= room_above ? below_y : above_y;
    }
    return clamp_position_to_bounds(
        QRect(QPoint(anchor_rect.left(), y), popup_size), bounds);
}

void show_popup_next_to_widget(QWidget* popup, QWidget* anchor)
{
    if (!popup->testAttribute(Qt::WA_Resized)) {
        // Match QWidget::show(), which sizes a window that was never explicitly
        // resized to fit its contents, so that we position the popup using the
        // size it will be shown at
        popup->adjustSize();
        popup->setAttribute(Qt::WA_Resized, false);
    }
    const QRect anchor_rect(anchor->mapToGlobal(QPoint(0, 0)), anchor->size());
    const auto bounds = anchor->screen()->availableGeometry();
    popup->move(calc_popup_position(popup->size(), anchor_rect, bounds));
    popup->show();
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
