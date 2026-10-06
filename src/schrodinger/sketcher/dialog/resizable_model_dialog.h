#pragma once

#include <vector>

#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace schrodinger
{
namespace sketcher
{

class DialogResizeHandle;

/**
 * Modal dialog with draggable edges and corners on WASM, where ModalDialog
 * removes the native window frame. Desktop builds use native window resizing.
 */
class SKETCHER_API ResizableModelDialog : public ModalDialog
{
    Q_OBJECT

  public:
    ResizableModelDialog(QWidget* parent = nullptr,
                         Qt::WindowFlags f = Qt::WindowFlags());

  protected:
    /**
     * Enable border handles for a frameless dialog. Called automatically on
     * WASM; native subclasses can also enable them when using a custom frame.
     */
    void enableBorderResizing();

    void resizeEvent(QResizeEvent* event) override;
    void showEvent(QShowEvent* event) override;

  private:
    void updateResizeHandles();

    // Owned by this dialog through QWidget parenting, outside its layout.
    std::vector<DialogResizeHandle*> m_resize_handles;
};

} // namespace sketcher
} // namespace schrodinger
