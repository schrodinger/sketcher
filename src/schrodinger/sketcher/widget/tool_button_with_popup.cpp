#include "schrodinger/sketcher/widget/tool_button_with_popup.h"

#include <stdexcept>

#include <QStyleOptionToolButton>
#include <QStylePainter>
#include <QTimer>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

ToolButtonWithPopup::ToolButtonWithPopup(QWidget* parent) : QToolButton(parent)
{
    setAttribute(Qt::WA_Hover, true);
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
    bool hovered = opt.state & QStyle::State_MouseOver;
    if (m_show_popup_indicator ||
        (m_show_popup_indicator_on_hover && hovered && isEnabled())) {
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
    configure_button_popup(this, popup_wdg);

    updateStyle();
}

void ToolButtonWithPopup::clearPopupWidget()
{
    m_popup_timer->stop();
    delete m_popup_wdg;
    m_popup_wdg = nullptr;
    showPopupIndicatorOnHover(false);
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
    m_show_popup_indicator = show;
    updateStyle();
}

void ToolButtonWithPopup::showPopupIndicatorOnHover(bool show)
{
    if (m_show_popup_indicator_on_hover == show) {
        return;
    }
    m_show_popup_indicator_on_hover = show;
    updateStyle();
    update();
}

void ToolButtonWithPopup::showPopup()
{
    if (m_popup_wdg == nullptr) {
        return;
    }

    show_popup_next_to_widget(m_popup_wdg, this);

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
    if (m_show_popup_indicator || m_show_popup_indicator_on_hover) {
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
