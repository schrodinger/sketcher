#include "schrodinger/sketcher/widget/widget_utils.h"

#include <algorithm>
#include <utility>

#include <QAbstractButton>
#include <QApplication>
#include <QButtonGroup>
#include <QEvent>
#include <QPointer>
#include <QScreen>

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

/**
 * @return The room available for a popup below and above an anchor
 */
static std::pair<int, int> calc_room_below_and_above(const QRect& anchor_rect,
                                                     const QRect& bounds)
{
    return {bounds.bottom() - anchor_rect.bottom(),
            anchor_rect.top() - bounds.top()};
}

int calc_max_popup_height(const QRect& anchor_rect, const QRect& bounds)
{
    const auto [room_below, room_above] =
        calc_room_below_and_above(anchor_rect, bounds);
    return std::max(room_below, room_above);
}

QPoint calc_popup_position(const QSize& popup_size, const QRect& anchor_rect,
                           const QRect& bounds)
{
    const auto [room_below, room_above] =
        calc_room_below_and_above(anchor_rect, bounds);
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
    popup->move(calc_popup_position(popup->size(), get_global_geometry(anchor),
                                    get_popup_bounds(anchor)));
    popup->show();
}

QRect get_global_geometry(const QWidget* widget)
{
    return QRect(widget->mapToGlobal(QPoint(0, 0)), widget->size());
}

QRect get_popup_bounds(const QWidget* anchor)
{
    return anchor->screen()->availableGeometry();
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
