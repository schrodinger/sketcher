#include "schrodinger/sketcher/dialog/resizable_model_dialog.h"

#include <algorithm>

#include <QMouseEvent>
#include <QResizeEvent>
#include <QShowEvent>

namespace schrodinger
{
namespace sketcher
{

static constexpr int RESIZE_HANDLE_WIDTH = 5;

/**
 * Invisible border widget that retains mouse input during a resize drag.
 */
class DialogResizeHandle : public QWidget
{
  public:
    DialogResizeHandle(Qt::Edges edges, QWidget* dialog) :
        QWidget(dialog),
        m_edges(edges)
    {
        setObjectName("dialog_resize_handle_" +
                      QString::number(static_cast<int>(edges)));
        const bool horizontal = edges & (Qt::LeftEdge | Qt::RightEdge);
        const bool vertical = edges & (Qt::TopEdge | Qt::BottomEdge);
        if (horizontal && vertical) {
            const bool forward =
                edges.testFlag(Qt::LeftEdge) == edges.testFlag(Qt::TopEdge);
            setCursor(forward ? Qt::SizeFDiagCursor : Qt::SizeBDiagCursor);
        } else {
            setCursor(horizontal ? Qt::SizeHorCursor : Qt::SizeVerCursor);
        }
    }

    void updatePosition()
    {
        const auto size = parentWidget()->size();
        const int border = RESIZE_HANDLE_WIDTH;
        int x = border;
        int y = border;
        int width = std::max(0, size.width() - 2 * border);
        int height = std::max(0, size.height() - 2 * border);
        if (m_edges & (Qt::LeftEdge | Qt::RightEdge)) {
            x = m_edges.testFlag(Qt::LeftEdge) ? 0 : size.width() - border;
            width = border;
        }
        if (m_edges & (Qt::TopEdge | Qt::BottomEdge)) {
            y = m_edges.testFlag(Qt::TopEdge) ? 0 : size.height() - border;
            height = border;
        }
        setGeometry(x, y, width, height);
        raise();
    }

  protected:
    void mousePressEvent(QMouseEvent* event) override
    {
        if (event->button() != Qt::LeftButton) {
            event->ignore();
            return;
        }
        m_dragging = true;
        m_start_position = event->globalPosition().toPoint();
        m_start_geometry = parentWidget()->geometry();
        event->accept();
    }

    void mouseMoveEvent(QMouseEvent* event) override
    {
        if (!m_dragging || !(event->buttons() & Qt::LeftButton)) {
            event->ignore();
            return;
        }
        auto* dialog = parentWidget();
        const auto delta = event->globalPosition().toPoint() - m_start_position;
        const auto hint = dialog->minimumSizeHint();
        const auto maximum = dialog->maximumSize();
        // An explicit minimum overrides the layout's minimum size hint,
        // matching QWidget's normal resizing rules.
        const auto minimum =
            QSize(dialog->minimumWidth() > 0 ? dialog->minimumWidth()
                                             : std::max(0, hint.width()),
                  dialog->minimumHeight() > 0 ? dialog->minimumHeight()
                                              : std::max(0, hint.height()))
                .boundedTo(maximum);
        auto geometry = m_start_geometry;
        if (m_edges & (Qt::LeftEdge | Qt::RightEdge)) {
            const bool left = m_edges.testFlag(Qt::LeftEdge);
            const int width = std::clamp(m_start_geometry.width() +
                                             (left ? -delta.x() : delta.x()),
                                         minimum.width(), maximum.width());
            geometry.setWidth(width);
            if (left) {
                geometry.moveLeft(m_start_geometry.right() - width + 1);
            }
        }
        if (m_edges & (Qt::TopEdge | Qt::BottomEdge)) {
            const bool top = m_edges.testFlag(Qt::TopEdge);
            const int height = std::clamp(m_start_geometry.height() +
                                              (top ? -delta.y() : delta.y()),
                                          minimum.height(), maximum.height());
            geometry.setHeight(height);
            if (top) {
                geometry.moveTop(m_start_geometry.bottom() - height + 1);
            }
        }
        dialog->setGeometry(geometry);
        event->accept();
    }

    void mouseReleaseEvent(QMouseEvent* event) override
    {
        if (event->button() == Qt::LeftButton) {
            m_dragging = false;
            event->accept();
        } else {
            event->ignore();
        }
    }

  private:
    Qt::Edges m_edges;
    bool m_dragging = false;
    QPoint m_start_position;
    QRect m_start_geometry;
};

ResizableModelDialog::ResizableModelDialog(QWidget* parent, Qt::WindowFlags f) :
    ModalDialog(parent, f)
{
#ifdef __EMSCRIPTEN__
    enableBorderResizing();
#endif
}

void ResizableModelDialog::enableBorderResizing()
{
    if (!m_resize_handles.empty()) {
        return;
    }
    const Qt::Edges edges[] = {Qt::LeftEdge,
                               Qt::RightEdge,
                               Qt::TopEdge,
                               Qt::BottomEdge,
                               Qt::LeftEdge | Qt::TopEdge,
                               Qt::RightEdge | Qt::TopEdge,
                               Qt::LeftEdge | Qt::BottomEdge,
                               Qt::RightEdge | Qt::BottomEdge};
    for (const auto edge : edges) {
        auto* handle = new DialogResizeHandle(edge, this);
        m_resize_handles.push_back(handle);
        handle->show();
    }
    updateResizeHandles();
}

void ResizableModelDialog::updateResizeHandles()
{
    for (auto* handle : m_resize_handles) {
        handle->updatePosition();
    }
}

void ResizableModelDialog::resizeEvent(QResizeEvent* event)
{
    ModalDialog::resizeEvent(event);
    updateResizeHandles();
}

void ResizableModelDialog::showEvent(QShowEvent* event)
{
    ModalDialog::showEvent(event);
    updateResizeHandles();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/resizable_model_dialog.moc"
