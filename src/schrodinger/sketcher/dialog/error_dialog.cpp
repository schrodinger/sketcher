#include "schrodinger/sketcher/dialog/error_dialog.h"

#include <QStyle>

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
    m_ui->text_edit->setText(text);
    m_ui->text_edit->setStyleSheet("QTextEdit { background: transparent; }");

    auto my_style = style();
    auto warning_icon = my_style->standardIcon(QStyle::SP_MessageBoxWarning);
    int icon_size = my_style->pixelMetric(QStyle::PM_MessageBoxIconSize);
    auto warning_pixmap = warning_icon.pixmap(icon_size, icon_size);
    m_ui->icon_lbl->setPixmap(warning_pixmap);
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