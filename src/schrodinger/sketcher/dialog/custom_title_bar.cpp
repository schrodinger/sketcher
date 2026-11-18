#include "schrodinger/sketcher/dialog/custom_title_bar.h"

#include <optional>

#include <QApplication>
#include <QColor>
#include <QHBoxLayout>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QPen>
#include <QToolButton>
#include <QtGlobal>

#include "schrodinger/sketcher/sketcher_css_style.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * The sizes, colors, and behavior of a CustomTitleBar
 */
struct TitleBarConfig {
    int height;
    int title_font_size;
    QColor background_color;
    // the color of the title bar's outline. If set to an invalid color (i.e.
    // QColor()), then no outline is painted
    QColor outline_color;
    // the font size of the close button, or std::nullopt if the title bar
    // shouldn't have a close button
    std::optional<int> close_button_font_size;
    // whether dragging the title bar moves the parent widget
    bool draggable;
};

static const TitleBarConfig DIALOG_TITLE_BAR_CONFIG(
    /* height: */ 30,
    /* title font size: */ 14,
    /* background color: */ QColor(64, 64, 64),
    /* outline color: */ QColor(Qt::black),
    // no close button
    /* close button font size: */ std::nullopt,
    /* draggable: */ true);

static const TitleBarConfig POPUP_TITLE_BAR_CONFIG(
    /* height: */ 24,
    /* title font size: */ 12,
    /* background color: */ QColor(POPUP_TITLE_BAR_COLOR),
    // no outline around the title bar
    /* outline color: */ QColor(),
    /* close button font size: */ 14,
    // the popup stays next to the widget that opened it
    /* draggable: */ false);

} // namespace

/**
 * @return The sizes, colors, and behavior for the given title bar type
 */
static const TitleBarConfig& get_config(TitleBarType type)
{
    if (type == TitleBarType::POPUP) {
        return POPUP_TITLE_BAR_CONFIG;
    } else if (type == TitleBarType::DIALOG) {
        return DIALOG_TITLE_BAR_CONFIG;
    }
    Q_UNREACHABLE_RETURN(DIALOG_TITLE_BAR_CONFIG);
}

CustomTitleBar::CustomTitleBar(const QString& title, QWidget* parent,
                               TitleBarType type) :
    QFrame(parent),
    m_type(type)
{
    const auto& config = get_config(type);
    QHBoxLayout* layout = new QHBoxLayout(this);
    this->setFixedHeight(config.height);
    m_title_label = new QLabel(title, this);
    m_title_label->setStyleSheet(
        QString("color: white; font-size: %1px; font-weight: 600;")
            .arg(config.title_font_size));
    layout->addWidget(m_title_label);
    if (config.close_button_font_size.has_value()) {
        addCloseButton(*config.close_button_font_size);
    }
}

CustomTitleBar::~CustomTitleBar() = default;

void CustomTitleBar::addCloseButton(int font_size)
{
    m_close_button = new QToolButton(this);
    m_close_button->setToolTip("Close");
    // a multiplication sign, as in the mockups
    m_close_button->setText(QString(QChar(0x00D7)));
    m_close_button->setStyleSheet(
        QString("QToolButton { color: white; font-size: %1px; "
                "    font-weight: bold; background-color: transparent; "
                "    border: none; }"
                "QToolButton:hover { color: #cccccc; }")
            .arg(font_size));
    auto* layout = qobject_cast<QHBoxLayout*>(this->layout());
    // Without vertical margins, the button fits in the fixed height. The title
    // is still centered vertically.
    auto margins = layout->contentsMargins();
    margins.setTop(0);
    margins.setBottom(0);
    margins.setRight(4);
    layout->setContentsMargins(margins);
    layout->addStretch();
    layout->addWidget(m_close_button);
}

QToolButton* CustomTitleBar::getCloseButton() const
{
    return m_close_button;
}

void CustomTitleBar::paintEvent(QPaintEvent* event)
{
    QPainter painter(this);
    const auto& config = get_config(m_type);
    painter.setBrush(config.background_color);
    QRectF paint_rect(rect());
    if (config.outline_color.isValid()) {
        const qreal pen_width = 1.0;
        painter.setPen(QPen(config.outline_color, pen_width));
        // A stroke is centered on the rectangle's edges, so we inset the
        // rectangle by half the pen width to keep the entire stroke inside of
        // the widget. Otherwise, part of the outline would be clipped, and
        // which edges are visible would depend on the device pixel ratio.
        const qreal inset = pen_width / 2;
        paint_rect.adjust(inset, inset, -inset, -inset);
    } else {
        painter.setPen(Qt::NoPen);
    }
    painter.drawRect(paint_rect);
}

void CustomTitleBar::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        m_cursor_start_pos = event->position();
        m_drag_started = false;
    }
}

void CustomTitleBar::mouseMoveEvent(QMouseEvent* event)
{
    if (!get_config(m_type).draggable) {
        return;
    }
    if (event->buttons() & Qt::LeftButton) {
        if (!m_drag_started &&
            (event->position() - m_cursor_start_pos).manhattanLength() >=
                QApplication::startDragDistance()) {
            m_drag_started = true;
        }
        if (m_drag_started && parentWidget()) {
            parentWidget()->move(parentWidget()->pos() +
                                 event->position().toPoint() -
                                 m_cursor_start_pos.toPoint());
        }
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/custom_title_bar.moc"
