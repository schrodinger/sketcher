#pragma once

#include <memory> // std::unique_ptr<>, required on Linux with Qt5

#include <QStyle>

#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace Ui
{
class MessageBoxDialog;
}

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API MessageBoxDialog : public ModalDialog
{
    Q_OBJECT
  public:
    MessageBoxDialog(const QString& title, const QString& text,
                     QStyle::StandardPixmap standard_icon,
                     QWidget* parent = nullptr,
                     Qt::WindowFlags f = Qt::WindowFlags());
    ~MessageBoxDialog();

  private:
    std::unique_ptr<Ui::MessageBoxDialog> m_ui;
};

/**
 * Convenience method for showing an error dialog
 */
void show_error_dialog(const QString& title, const QString& text,
                       QWidget* parent, Qt::WindowFlags f = Qt::WindowFlags());

/**
 * Convenience method for showing an information dialog
 */
void show_information_dialog(const QString& title, const QString& text,
                             QWidget* parent,
                             Qt::WindowFlags f = Qt::WindowFlags());

} // namespace sketcher
} // namespace schrodinger
