#include "schrodinger/sketcher/dialog/modal_dialog.h"

#include <QVBoxLayout>
#include <QWidget>

namespace schrodinger
{
namespace sketcher
{

ModalDialog::ModalDialog(QWidget* parent, Qt::WindowFlags f) :
    QDialog(parent, f)
{
    setWindowModality(Qt::ApplicationModal);
    setAttribute(Qt::WA_DeleteOnClose, true);

    m_dlg_layout = new QVBoxLayout();

// Add a custom title bar to the dialog if we're in WASM
#ifdef __EMSCRIPTEN__
    // Remove the default title bar drawn by emscripten
    // Note: this has the side effect of removing the window's border
    setWindowFlags(Qt::FramelessWindowHint | Qt::Dialog);

    // Create a wrapper layout that holds the title bar and the dialog contents
    // This layout should have 0 margins for styling purposes
    QVBoxLayout* wrapper_layout = new QVBoxLayout(this);
    wrapper_layout->addLayout(m_dlg_layout);
    this->setLayout(wrapper_layout);
    wrapper_layout->setContentsMargins(0, 0, 0, 0);
    // Removing the window's border makes the dialog contents look closer
    // to the edge of the window, so we add some margins to the dialog contents
    // to compensate
    m_dlg_layout->setContentsMargins(5, 5, 5, 5);

    m_title_bar = new CustomTitleBar("Placeholder Title", this);
    wrapper_layout->insertWidget(0, m_title_bar);

    // When the window title changes, update the title bar's title
    connect(this, &QDialog::windowTitleChanged, this,
            &ModalDialog::onWindowTitleChanged);
#endif
}

ModalDialog::~ModalDialog() = default;

void ModalDialog::onWindowTitleChanged(const QString& newTitle)
{
    if (m_title_bar != nullptr) {
        m_title_bar->m_title_label->setText(newTitle);
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/modal_dialog.moc"