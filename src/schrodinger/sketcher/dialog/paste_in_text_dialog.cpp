#include "paste_in_text_dialog.h"
#include "schrodinger/sketcher/ui/ui_paste_in_text_dialog.h"
#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

PasteInTextDialog::PasteInTextDialog(SketcherModel* model, QWidget* parent) :
    ModalDialog(parent)
{
    if (model == nullptr) {
        throw std::runtime_error(
            "This dialog cannot be created without a model.");
    }
    m_sketcher_model = model;

    ui.reset(new Ui::PasteInTextDialog());
    ui->setupUi(this);

    connect(model, &SketcherModel::valuesChanged, this,
            &PasteInTextDialog::onModelValuesChanged);
    onModelValuesChanged();
}

PasteInTextDialog::~PasteInTextDialog() = default;

void PasteInTextDialog::onModelValuesChanged()
{
    bool replace =
        m_sketcher_model->getValue(ModelKey::NEW_STRUCTURES_REPLACE_CONTENT)
            .toBool();
    QString text =
        replace ? "Specified structure will <b>replace</b> Sketcher content"
                : "Specified structure will be added to Sketcher content";
    QString style = replace ? "color : #c87c00" : "color : gray";
    ui->status_label->setText(text);
    ui->status_label->setStyleSheet(style);
}

void PasteInTextDialog::accept()
{
    QDialog::accept();
    emit textAccepted(ui->structure_text_edit->toPlainText());
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/paste_in_text_dialog.moc"
