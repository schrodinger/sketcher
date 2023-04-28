#pragma once

#include <QFrame>
#include <QLabel>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Custom title bar with a fixed height and styled QLabel for the title.
 * The default title bar drawn by emscripten looks really bad so we draw our
 * own. We override the mouse events to allow the user to drag the parent
 * widget. Note: For styling purposes, parent widget's layout needs to have a
 * margin of 0.
 */
class SKETCHER_API CustomTitleBar : public QFrame
{
    Q_OBJECT
  public:
    CustomTitleBar(const QString& title, QWidget* parent = nullptr);
    ~CustomTitleBar();
    QLabel* m_title_label = nullptr;

  private:
    /**
     * Paints the background of the title bar.
     */
    void paintEvent(QPaintEvent* event) override;

    /**
     * Sets the start position of the mouse press event.
     */
    void mousePressEvent(QMouseEvent* event) override;

    /**
     * Moves the parent widget to the new position.
     */
    void mouseMoveEvent(QMouseEvent* event) override;

    QPointF m_cursor_start_pos;
    bool m_drag_started;
};

} // namespace sketcher
} // namespace schrodinger
