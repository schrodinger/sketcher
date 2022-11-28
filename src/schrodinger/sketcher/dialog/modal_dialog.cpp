#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace schrodinger
{
namespace sketcher
{

ModalDialog::ModalDialog(QWidget* parent, Qt::WindowFlags f) :
    QDialog(parent, f)
{
    setWindowModality(Qt::ApplicationModal);
    setAttribute(Qt::WA_DeleteOnClose, true);
}

ModalDialog::~ModalDialog() = default;

void ModalDialog::paintEvent(QPaintEvent* event)
{
    // SKETCH-1600: Alter dialog appearance if used in the browser
    QDialog::paintEvent(event);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/modal_dialog.moc"