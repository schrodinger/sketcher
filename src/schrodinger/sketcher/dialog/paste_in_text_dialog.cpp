#include "paste_in_text_dialog.h"

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_paste_in_text_dialog.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{
namespace
{

/**
 * @return a list of tuples of
 *   - format name
 *   - format enum
 *   - whether the format applies to atomistic models, monomeric models, or both
 * for all formats to load into the format combo box
 */
static std::vector<std::tuple<QString, Format, InterfaceTypeType>>
build_format_list()
{
    std::vector<std::tuple<QString, Format, InterfaceTypeType>> formats = {
        {"Autodetect", Format::AUTO_DETECT,
         InterfaceType::ATOMISTIC_OR_MONOMERIC},
    };
    auto append_to_formats =
        [&formats](std::vector<std::tuple<Format, std::string>> to_append,
                   InterfaceTypeType fmt_mol_type) {
            std::transform(
                to_append.begin(), to_append.end(), std::back_inserter(formats),
                [fmt_mol_type](auto cur_format)
                    -> std::tuple<QString, Format, InterfaceTypeType> {
                    auto [fmt_val, fmt_name] = cur_format;
                    return {QString::fromStdString(fmt_name), fmt_val,
                            fmt_mol_type};
                });
        };
    auto mol_import_formats = get_mol_import_formats();
    auto monomeric_import_formats = get_monomoric_import_formats();
    append_to_formats(mol_import_formats, InterfaceType::ATOMISTIC);
    append_to_formats(monomeric_import_formats, InterfaceType::MONOMERIC);
    return formats;
}

/**
 * All formats to load into the format combo box, and whether they apply to
 * atomistic or monomeric models
 */
const auto COMBO_BOX_FORMATS = build_format_list();

} // namespace

PasteInTextDialog::PasteInTextDialog(SketcherModel* model, QWidget* parent) :
    ModalDialog(parent)
{
    if (model == nullptr) {
        throw std::runtime_error(
            "This dialog cannot be created without a model.");
    }
    m_sketcher_model = model;

    ui.reset(new Ui::PasteInTextDialog());
    setupDialogUI(*ui);

    // Connect signals and slots
    connect(model, &SketcherModel::valuesChanged, this,
            &PasteInTextDialog::onModelValuesChanged);
    connect(ui->structure_text_edit, &QPlainTextEdit::textChanged, this,
            &PasteInTextDialog::autodetectCurrentFormat);
    onModelValuesChanged();
}

PasteInTextDialog::~PasteInTextDialog() = default;

void PasteInTextDialog::onModelValuesChanged()
{
    bool replace = m_sketcher_model->getValueBool(
        ModelKey::NEW_STRUCTURES_REPLACE_CONTENT);
    QString text =
        replace ? "Specified structure will <b>replace</b> Sketcher content"
                : "Specified structure will be added to Sketcher content";
    QString style = replace ? "color : #c87c00" : "color : gray";
    ui->status_lbl->setText(text);
    ui->status_lbl->setStyleSheet(style);

    // load data into the format combo box
    auto interface_type = m_sketcher_model->getInterfaceType();
    auto cur_mol_type = m_sketcher_model->getMoleculeType();
    // first figure out whether we want atomistic formats, monomeric formats, or
    // both
    auto allowed_mol_type = InterfaceType::ATOMISTIC_OR_MONOMERIC;
    if (interface_type == InterfaceType::ATOMISTIC ||
        (!replace && cur_mol_type == MoleculeType::ATOMISTIC)) {
        allowed_mol_type = InterfaceType::ATOMISTIC;
    } else if (interface_type == InterfaceType::MONOMERIC ||
               (!replace && cur_mol_type == MoleculeType::MONOMERIC)) {
        allowed_mol_type = InterfaceType::MONOMERIC;
    }
    // then iterate through the available formats and load tha applicable ones
    // into the combo box
    ui->format_combo->clear();
    for (const auto& [format_name, format, format_mol_type] :
         COMBO_BOX_FORMATS) {
        if (!(format_mol_type & allowed_mol_type)) {
            // this format doesn't apply, so skip it
            continue;
        }
        if (format == Format::AUTO_DETECT) {
            // remember where the autodetect item is so we can update the name
            // later with the autodetected format
            m_autodetect_index = ui->format_combo->count();
            m_autodetect_base_name = format_name;
        } else {
            m_formats[format] = format_name;
        }
        ui->format_combo->addItem(format_name, QVariant::fromValue(format));
    }
    autodetectCurrentFormat();
}

/**
 * Determine whether the given text can be parsed as the specified format
 */
static bool text_is_readable_as(const std::string& text, const Format format)
{
    try {
        convert_text_to_mol_or_reaction(text, format);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

void PasteInTextDialog::autodetectCurrentFormat()
{
    auto qtext = ui->structure_text_edit->toPlainText();
    auto text = qtext.toStdString();
    bool found_format = false;
    if (qtext.trimmed().isEmpty()) {
        m_autodetected_format = std::nullopt;
        found_format = true;
    } else {
        for (const auto cur_format : rdkit_extensions::AUTO_DETECT_FORMATS) {
            if (!m_formats.contains(cur_format)) {
                continue;
            }
            if (text_is_readable_as(text, cur_format)) {
                m_autodetected_format = cur_format;
                found_format = true;
                break;
            }
        }
    }
    QString autodetect_label = m_autodetect_base_name;
    if (!found_format) {
        // the user has entered text, but we can't parse it
        m_autodetected_format = std::nullopt;
        autodetect_label += " (unrecognized)";
    } else if (m_autodetected_format.has_value()) {
        // there is text and we could parse it
        auto format_name = m_formats.at(*m_autodetected_format);
        autodetect_label += QString(" (%1)").arg(format_name);
    }
    ui->format_combo->setItemText(m_autodetect_index, autodetect_label);
}

void PasteInTextDialog::accept()
{
    QDialog::accept();
    auto format = ui->format_combo->currentData().value<Format>();
    if (format == Format::AUTO_DETECT && m_autodetected_format.has_value()) {
        format = *m_autodetected_format;
    }
    emit textAccepted(ui->structure_text_edit->toPlainText().toStdString(),
                      format);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/paste_in_text_dialog.moc"
