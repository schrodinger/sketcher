#include "schrodinger/sketcher/dialog/message_box_dialog.h"

#include "schrodinger/sketcher/ui/ui_message_box_dialog.h"

namespace schrodinger
{
namespace sketcher
{

MessageBoxDialog::MessageBoxDialog(const QString& title, const QString& text,
                                   QStyle::StandardPixmap standard_icon,
                                   QWidget* parent, Qt::WindowFlags f) :
    ModalDialog(parent, f)
{
    m_ui.reset(new Ui::MessageBoxDialog());
    setupDialogUI(*m_ui);

    setWindowTitle(title);
    m_ui->text_edit->setText(text);
    m_ui->text_edit->setStyleSheet("QTextEdit { background: transparent; }");

    auto my_style = style();
    auto icon = my_style->standardIcon(standard_icon);
    int icon_size = my_style->pixelMetric(QStyle::PM_MessageBoxIconSize);
    auto pixmap = icon.pixmap(icon_size, icon_size);
    m_ui->icon_lbl->setPixmap(pixmap);
}

MessageBoxDialog::~MessageBoxDialog() = default;

void show_error_dialog(const QString& title, const QString& text,
                       QWidget* parent, Qt::WindowFlags f)
{
    auto error_dlg = new MessageBoxDialog(
        title, text, QStyle::SP_MessageBoxWarning, parent, f);
    error_dlg->show();
}

void show_information_dialog(const QString& title, const QString& text,
                             QWidget* parent, Qt::WindowFlags f)
{
    auto information_dlg = new MessageBoxDialog(
        title, text, QStyle::SP_MessageBoxInformation, parent, f);
    information_dlg->show();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/message_box_dialog.moc"
