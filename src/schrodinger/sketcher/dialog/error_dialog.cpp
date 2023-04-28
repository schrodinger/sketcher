#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/ui/ui_error_dialog.h"

namespace schrodinger
{
namespace sketcher
{

ErrorDialog::ErrorDialog(const QString& title, const QString& text,
                         QWidget* parent, Qt::WindowFlags f) :
    ModalDialog(parent, f)
{
    m_ui.reset(new Ui::ErrorDialog());
    setupDialogUI(*m_ui);

    setWindowTitle(title);
    m_ui->text_lbl->setText(text);
}

ErrorDialog::~ErrorDialog() = default;

void show_error_dialog(const QString& title, const QString& text,
                       QWidget* parent, Qt::WindowFlags f)
{
    auto error_dlg = new ErrorDialog(title, text, parent, f);
    error_dlg->show();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/error_dialog.moc"