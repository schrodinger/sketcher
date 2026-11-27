#pragma once

#include <QFrame>
#include <QLabel>

#include "schrodinger/sketcher/definitions.h"

class QToolButton;

namespace schrodinger
{
namespace sketcher
{

/**
 * The kind of a CustomTitleBar, which determines how it's drawn and how it
 * behaves. Both types use white text on a dark background.
 */
enum class TitleBarType {
    // a dialog's title bar, which replaces emscripten's title bar on WASM. It
    // can be dragged to move the dialog and has no close button.
    DIALOG,
    // a smaller, lighter title bar for popups on any platform. It can't be
    // dragged, since popups stay next to the widget that opened them, and it
    // has a close button.
    POPUP,
};

/**
 * Custom title bar with a fixed height and styled QLabel for the title. The
 * default title bar drawn by emscripten looks really bad so dialogs draw their
 * own on WASM. Popups never have a native title bar, so popups with a title bar
 * use this on any platform.
 *
 * Note: For styling purposes, parent widget's layout needs to have a margin of
 * 0.
 */
class SKETCHER_API CustomTitleBar : public QFrame
{
    Q_OBJECT
  public:
    CustomTitleBar(const QString& title, QWidget* parent = nullptr,
                   TitleBarType type = TitleBarType::DIALOG);
    ~CustomTitleBar();
    QLabel* m_title_label = nullptr;

    /**
     * @return The close button, or nullptr if this type of title bar doesn't
     * have one. The title bar doesn't connect it to anything.
     */
    QToolButton* getCloseButton() const;

  private:
    /**
     * Add a close button to the right end of the title bar
     *
     * @param font_size The font size of the close button, in pixels
     */
    void addCloseButton(int font_size);

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

    TitleBarType m_type;
    QToolButton* m_close_button = nullptr;
    QPointF m_cursor_start_pos;
    bool m_drag_started = false;
};

} // namespace sketcher
} // namespace schrodinger
