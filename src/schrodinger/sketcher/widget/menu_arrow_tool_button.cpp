#include "schrodinger/sketcher/widget/menu_arrow_tool_button.h"

#include <QStyleOptionButton>
#include <QStyleOptionToolButton>
#include <QStylePainter>

namespace schrodinger
{
namespace sketcher
{

MenuArrowToolButton::MenuArrowToolButton(QWidget* parent) : QToolButton(parent)
{
    setAttribute(Qt::WA_Hover, true);
}

void MenuArrowToolButton::paintEvent(QPaintEvent* event)
{
    QStylePainter painter(this);
    QStyleOptionToolButton option;
    initStyleOption(&option);
    if (isChecked() ||
        (isEnabled() && (option.state & QStyle::State_MouseOver))) {
        option.features |= QStyleOptionToolButton::HasMenu;
    }

    // Paint the full-row background and native menu indicator separately from
    // the label so the icon and text stay centered together.
    auto background = option;
    background.text.clear();
    background.icon = QIcon();
    painter.drawComplexControl(QStyle::CC_ToolButton, background);

    QStyleOptionButton label;
    label.initFrom(this);
    label.rect = option.rect;
    label.state = option.state;
    label.palette = option.palette;
    label.fontMetrics = option.fontMetrics;
    label.text = option.text;
    label.icon = option.icon;
    label.iconSize = option.iconSize;
    painter.setFont(option.font);
    painter.drawControl(QStyle::CE_PushButtonLabel, label);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/menu_arrow_tool_button.moc"
