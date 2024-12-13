#include "schrodinger/sketcher/dialog/file_export_dialog.h"

#include <QFileDialog>
#include <QTimer>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_file_export_dialog.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

namespace
{

const QString DEFAULT_FILENAME{"structure"};

const auto get_format_list = [](bool has_reaction) {
    return has_reaction ? get_reaction_export_formats()
                        : get_standard_export_formats();
};

} // namespace

FileExportDialog::FileExportDialog(SketcherModel* model, QWidget* parent) :
    ModalDialog(parent),
    m_model_has_reaction(model->hasReaction())
{
    m_ui.reset(new Ui::FileExportDialog());
    setupDialogUI(*m_ui);

#ifdef __EMSCRIPTEN__
    // In the Browser the exporting of data is a Download action.
    // In Maestro the exporting of data is a File Save action.
    // The button should say "Download" with no "...".
    // In Maestro, it should continue to say "Save...".
    m_ui->export_btn->setText("Download");
#endif
    // Connect signals and slots
    connect(m_ui->export_btn, &QPushButton::clicked, this,
            &FileExportDialog::exportFile);
    connect(m_ui->cancel_btn, &QPushButton::clicked, this, &QDialog::reject);

    m_ui->filename_le->setText(DEFAULT_FILENAME);

    // initialize the format combo
    for (const auto& [fmt, label, extensions] :
         get_format_list(m_model_has_reaction)) {
        if (!extensions.empty()) {
            auto filter = get_filter_name(label, extensions);
            m_ui->format_combo->addItem(filter, QVariant::fromValue(fmt));
        }
    }
}

FileExportDialog::~FileExportDialog() = default;

QByteArray FileExportDialog::getFileContent() const
{
    auto combo_format = m_ui->format_combo->currentData().value<Format>();
    return emit exportTextRequested(combo_format).toUtf8();
}

std::vector<std::string> FileExportDialog::getValidExtensions() const
{
    auto combo_format = m_ui->format_combo->currentData().value<Format>();
    return m_model_has_reaction ? get_rxn_extensions(combo_format)
                                : get_mol_extensions(combo_format);
}

void FileExportDialog::exportFile()
{
    // Ensure filename hint has an appropriate extension
    auto filename = m_ui->filename_le->text();
    auto suffix = QFileInfo(filename).completeSuffix();
    auto extensions = getValidExtensions();
    auto suffix_not_found =
        std::find(extensions.begin(), extensions.end(),
                  "." + suffix.toStdString()) == extensions.end();
    if (suffix.isEmpty() || suffix_not_found) {
        filename.append(extensions.at(0));
    }

    QByteArray file_content;
    try {
        file_content = getFileContent();
    } catch (const std::exception& exc) {
        show_error_dialog("Export Error", exc.what(), parentWidget());
        reject();
        return;
    }

    using namespace rdkit_extensions;

    if (auto compression_type =
            get_compression_type_from_ext(filename.toStdString());
        compression_type != CompressionType::UNKNOWN) {
        auto text =
            get_compressed_string(file_content.toStdString(), compression_type);
        file_content = QByteArray::fromStdString(text);
    }

    QFileDialog::saveFileContent(file_content, filename);
    accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/file_export_dialog.moc"
