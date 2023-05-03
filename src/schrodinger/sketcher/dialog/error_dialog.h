#pragma once

#include <memory> // std::unique_ptr<>, required on Linux with Qt5

#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace Ui
{
class ErrorDialog;
}

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ErrorDialog : public ModalDialog
{
    Q_OBJECT
  public:
    ErrorDialog(const QString& title, const QString& text,
                QWidget* parent = nullptr,
                Qt::WindowFlags f = Qt::WindowFlags());
    ~ErrorDialog();

  private:
    std::unique_ptr<Ui::ErrorDialog> m_ui;
};

/**
 * Convience method for showing an error dialog
 */
void show_error_dialog(const QString& title, const QString& text,
                       QWidget* parent, Qt::WindowFlags f = Qt::WindowFlags());

} // namespace sketcher
} // namespace schrodinger
