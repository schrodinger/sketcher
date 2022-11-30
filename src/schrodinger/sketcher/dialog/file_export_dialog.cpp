#include "schrodinger/sketcher/dialog/file_export_dialog.h"

#include <QFileDialog>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/file_import_export.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_file_export_dialog.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

namespace
{

const QString DEFAULT_FILENAME{"structure"};

const auto get_format_list = [](bool has_reaction) {
    return has_reaction ? REACTION_FORMATS : STANDARD_FORMATS;
};

} // namespace

FileExportDialog::FileExportDialog(SketcherModel* model, QWidget* parent) :
    ModalDialog(parent),
    m_model_has_reaction(model->hasReaction())
{
    m_ui.reset(new Ui::FileExportDialog());
    m_ui->setupUi(this);

    connect(m_ui->export_btn, &QPushButton::clicked, this,
            &FileExportDialog::exportFile);

    m_ui->filename_le->setText(DEFAULT_FILENAME);
    // initialize the format combo
    auto format_list = get_format_list(m_model_has_reaction);
    for (const auto& [format, filter] : get_name_filters(format_list)) {
        m_ui->format_combo->addItem(filter, QVariant::fromValue(format));
    }
}

FileExportDialog::~FileExportDialog() = default;

QByteArray FileExportDialog::getFileContent() const
{
    auto combo_format = m_ui->format_combo->currentData().value<Format>();
    return emit exportTextRequested(combo_format).toUtf8();
}

QStringList FileExportDialog::getValidExtensions() const
{
    auto format_list = get_format_list(m_model_has_reaction);
    auto combo_format = m_ui->format_combo->currentData().value<Format>();
    return get_file_extensions(format_list, combo_format);
}

void FileExportDialog::exportFile()
{
    // Ensure filename hint has an appropriate extension
    auto filename = m_ui->filename_le->text();
    auto suffix = QFileInfo(filename).completeSuffix();
    auto extensions = getValidExtensions();
    if (suffix.isEmpty() || !extensions.contains("." + suffix)) {
        filename.append(extensions.at(0));
    }

    QByteArray file_content;
    try {
        file_content = getFileContent();
    } catch (const std::exception& exc) {
        auto text = QString("Cannot export: ") + exc.what();
        show_error_dialog("Export Error", text, parentWidget());
        reject();
        return;
    }

    QFileDialog::saveFileContent(file_content, filename);
    accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/file_export_dialog.moc"
