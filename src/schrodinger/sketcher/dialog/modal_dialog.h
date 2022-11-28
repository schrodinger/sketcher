#pragma once

#include <QDialog>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ModalDialog : public QDialog
{
    Q_OBJECT
  public:
    ModalDialog(QWidget* parent = nullptr,
                Qt::WindowFlags f = Qt::WindowFlags());
    ~ModalDialog();

  private:
    void paintEvent(QPaintEvent* event) override;

    // NOTE: modal dialogs that use `exec()` are problematic in the
    // WebAssembly (WASM) build, and will not return values. Explicitly enforce
    // dialogs inherited from this class to be modal and use show.
    using QDialog::exec;
    using QDialog::open;
};

} // namespace sketcher
} // namespace schrodinger
