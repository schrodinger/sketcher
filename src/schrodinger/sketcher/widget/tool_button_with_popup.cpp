#include "schrodinger/sketcher/widget/tool_button_with_popup.h"

#include <stdexcept>

#include <QApplication>
#include <QPoint>
#include <QScreen>
#include <QStyleOptionToolButton>
#include <QStylePainter>
#include <QTimer>

#include "schrodinger/sketcher/sketcher_css_style.h"

namespace schrodinger
{
namespace sketcher
{

ToolButtonWithPopup::ToolButtonWithPopup(QWidget* parent) : QToolButton(parent)
{
    m_popup_timer = new QTimer(this);
    m_popup_timer->setSingleShot(true);
    m_popup_timer->setInterval(m_popup_delay);
    connect(m_popup_timer, &QTimer::timeout, this,
            &ToolButtonWithPopup::showPopup);
    connect(this, &ToolButtonWithPopup::pressed, m_popup_timer,
            static_cast<void (QTimer::*)()>(&QTimer::start));
    connect(this, &ToolButtonWithPopup::released, m_popup_timer, &QTimer::stop);
    connect(this, &ToolButtonWithPopup::pressed, this,
            &ToolButtonWithPopup::onPressed);
    connect(this, &ToolButtonWithPopup::clicked, this,
            &ToolButtonWithPopup::onClicked);
}

void ToolButtonWithPopup::paintEvent(QPaintEvent* event)
{
    QStylePainter p(this);
    QStyleOptionToolButton opt;
    initStyleOption(&opt);
    if (m_show_corner_arrow) {
        opt.features |= QStyleOptionToolButton::HasMenu;
    }
    p.drawComplexControl(QStyle::CC_ToolButton, opt);
}

void ToolButtonWithPopup::setPopupWidget(QWidget* popup_wdg)
{
    if (m_popup_wdg != nullptr) {
        throw std::runtime_error("The popup widget has already been set.");
    }
    m_popup_wdg = popup_wdg;

    updateStyle();
}

QWidget* ToolButtonWithPopup::getPopupWidget() const
{
    return m_popup_wdg;
}

void ToolButtonWithPopup::setPopupDelay(float popup_delay)
{
    if (m_popup_delay != popup_delay) {
        m_popup_delay = popup_delay;
        m_popup_timer->setInterval(m_popup_delay);
    }
}

void ToolButtonWithPopup::showPopupIndicator(bool show)
{
    m_show_corner_arrow = show;
    updateStyle();
}

void ToolButtonWithPopup::showPopup()
{
    if (m_popup_wdg == nullptr) {
        return;
    }

    // The coordinates at the top left of this button
    auto button_pos = mapToGlobal(QPoint());
    auto button_x = button_pos.x();
    auto button_y = button_pos.y();
    // A rectangle for the desktop screen
    auto screen_rec = QApplication::screenAt(QCursor::pos())->geometry();
    auto screen_x = screen_rec.x();
    auto screen_y = screen_rec.y();
    // A rectangle for the popup
    auto popup_rec = m_popup_wdg->geometry();

    float pos_x = 0;
    if (button_x < screen_x) {
        // Edge case where cursor is not on the same screen as the top left
        // corner of the button
        pos_x = screen_x;
    } else {
        pos_x = button_x + std::min(0, screen_x + screen_rec.width() -
                                           button_x - popup_rec.width());
    }

    // Include the height of the button because we want to show the popup below
    // the button by default
    float pos_y = 0;
    if (button_y < screen_y) {
        // Edge case where cursor is not on the same screen as the top left
        // corner of the button
        pos_y = screen_y;
    } else if (button_y + popup_rec.height() + geometry().height() >
               screen_y + screen_rec.height()) {
        pos_y = button_y - popup_rec.height();
    } else {
        pos_y = button_y + geometry().height();
    }
    m_popup_wdg->move(QPoint(pos_x, pos_y));
    m_popup_wdg->show();

    // If the user presses and holds to show the popup, we want this button to
    // be activated as if it were clicked
    if (isCheckable() && !isChecked()) {
        click();
    }
}

void ToolButtonWithPopup::setStyleSheet(QString text)
{
    m_custom_style_sheet = text;
    updateStyle();
}

void ToolButtonWithPopup::updateStyle()
{
    QString style;
    if (m_show_corner_arrow) {
        style.append(TOOL_BUTTON_CORNER_ARROW_STYLE);
    }
    style.append(m_custom_style_sheet);
    QToolButton::setStyleSheet(style);
}

void ToolButtonWithPopup::onPressed()
{
    m_was_checked_on_press = isCheckable() && isChecked();
}

void ToolButtonWithPopup::onClicked(bool checked)
{
    if (m_was_checked_on_press && m_popup_wdg != nullptr &&
        !m_popup_wdg->isVisibleTo(this)) {
        // If the button was checked when the user clicked it, show popup
        showPopup();
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/tool_button_with_popup.moc"
