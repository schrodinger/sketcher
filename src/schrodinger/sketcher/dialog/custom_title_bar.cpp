#include "schrodinger/sketcher/dialog/custom_title_bar.h"

#include <QApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>

namespace schrodinger
{
namespace sketcher
{

CustomTitleBar::CustomTitleBar(const QString& title, QWidget* parent) :
    QFrame(parent)
{
    QHBoxLayout* layout = new QHBoxLayout(this);
    this->setFixedHeight(30);
    m_title_label = new QLabel(title, this);
    m_title_label->setStyleSheet(
        "color: white; font-size: 14px; font-weight: 600;");
    layout->addWidget(m_title_label);
}

CustomTitleBar::~CustomTitleBar() = default;

void CustomTitleBar::paintEvent(QPaintEvent* event)
{
    QPainter painter(this);
    painter.setBrush(QBrush(QColor(64, 64, 64)));
    painter.drawRect(rect());
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