#include "schrodinger/sketcher/dialog/modal_dialog.h"

#include <QCheckBox>
#include <QComboBox>
#include <QLineEdit>
#include <QRadioButton>
#include <QSpinBox>
#include <QTimer>
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

    // Create a wrapper layout that holds the title bar and the dialog contents
    QVBoxLayout* wrapper_layout = new QVBoxLayout(this);
    setLayout(wrapper_layout);
    m_dlg_layout = new QVBoxLayout();
    wrapper_layout->addLayout(m_dlg_layout);

// Add a custom title bar to the dialog if we're in WASM
#ifdef __EMSCRIPTEN__
    // Remove the default title bar drawn by emscripten
    // Note: this has the side effect of removing the window's border
    setWindowFlags(Qt::FramelessWindowHint | Qt::Dialog);

    // This layout should have 0 margins for styling purposes
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

void connect_input_widgets_to_timer(const QWidget* const dialog,
                                    const QTimer* const timer)
{
    QList<QObject*> descendents = dialog->children();
    while (!descendents.isEmpty()) {
        auto* child = descendents.takeFirst();
        if (auto* combo = qobject_cast<QComboBox*>(child)) {
            dialog->connect(combo, &QComboBox::currentIndexChanged, timer,
                            qOverload<>(&QTimer::start));
        } else if (auto* sb = qobject_cast<QSpinBox*>(child)) {
            dialog->connect(sb, &QSpinBox::valueChanged, timer,
                            qOverload<>(&QTimer::start));
        } else if (auto* sb = qobject_cast<QDoubleSpinBox*>(child)) {
            dialog->connect(sb, &QDoubleSpinBox::valueChanged, timer,
                            qOverload<>(&QTimer::start));
        } else if (auto* rb = qobject_cast<QRadioButton*>(child)) {
            dialog->connect(rb, &QRadioButton::toggled, timer,
                            qOverload<>(&QTimer::start));
        } else if (auto* cb = qobject_cast<QCheckBox*>(child)) {
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
            dialog->connect(cb, &QCheckBox::checkStateChanged, timer,
                            qOverload<>(&QTimer::start));
#else
            dialog->connect(cb, &QCheckBox::stateChanged, timer,
                            qOverload<>(&QTimer::start));
#endif
        } else if (auto* le = qobject_cast<QLineEdit*>(child)) {
            dialog->connect(le, &QLineEdit::textChanged, timer,
                            qOverload<>(&QTimer::start));
        } else {
            descendents.append(child->children());
        }
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/modal_dialog.moc"
