#include "paste_in_text_dialog.h"

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_paste_in_text_dialog.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{

namespace
{

typedef uint8_t FormatMoleculeType_t;
namespace FormatMoleculeType
{
/**
 * Whether a given format applies to atomistic models, monomeric models, or both
 */
enum : FormatMoleculeType_t {
    ATOMISTIC = 1 << 0,
    MONOMERIC = 1 << 1,
    BOTH = ATOMISTIC | MONOMERIC,
};
} // namespace FormatMoleculeType

/**
 * All formats to load into the format combo box, and whether they apply to
 * atomistic or monomeric models
 */
const std::vector<std::tuple<QString, Format, FormatMoleculeType_t>> FORMATS = {
    {"Autodetect", Format::AUTO_DETECT, FormatMoleculeType::BOTH},
    {"RDKit Base64", Format::RDMOL_BINARY_BASE64, FormatMoleculeType::BOTH},
    {"SMILES", Format::SMILES, FormatMoleculeType::ATOMISTIC},
    {"Extended SMILES", Format::EXTENDED_SMILES, FormatMoleculeType::ATOMISTIC},
    {"SMARTS", Format::SMARTS, FormatMoleculeType::ATOMISTIC},
    {"Extended SMARTS", Format::EXTENDED_SMARTS, FormatMoleculeType::ATOMISTIC},
    {"MDL SD V3000", Format::MDL_MOLV3000, FormatMoleculeType::ATOMISTIC},
    {"MDL SD V2000", Format::MDL_MOLV2000, FormatMoleculeType::ATOMISTIC},
    {"Maestro", Format::MAESTRO, FormatMoleculeType::ATOMISTIC},
    {"InChI", Format::INCHI, FormatMoleculeType::ATOMISTIC},
    {"PDB", Format::PDB, FormatMoleculeType::ATOMISTIC},
    {"MOL2", Format::MOL2, FormatMoleculeType::ATOMISTIC},
    {"XYZ", Format::XYZ, FormatMoleculeType::ATOMISTIC},
    {"Marvin Document", Format::MRV, FormatMoleculeType::ATOMISTIC},
    {"ChemDraw XML", Format::CDXML, FormatMoleculeType::ATOMISTIC},
    {"HELM", Format::HELM, FormatMoleculeType::MONOMERIC},
    {"FASTA Peptide", Format::FASTA_PEPTIDE, FormatMoleculeType::MONOMERIC},
    {"FASTA DNA", Format::FASTA_DNA, FormatMoleculeType::MONOMERIC},
    {"FASTA RNA", Format::FASTA_RNA, FormatMoleculeType::MONOMERIC},
};

/**
 * Formats to auto-detect, in the order that they'll be tried. Note that we'll
 * skip any atomistic/monomeric formats here if the SketcherWidget can't
 * currently accept that type of model.
 */
const std::vector<Format> AUTO_DETECT_FORMATS = {
    Format::RDMOL_BINARY_BASE64,
    Format::MDL_MOLV3000,
    Format::MDL_MOLV2000,
    Format::MAESTRO,
    Format::INCHI,
    Format::PDB,
    Format::MOL2,
    Format::XYZ,
    Format::MRV,
    Format::CDXML,
    // Attempt SMILES before SMARTS, given not all SMARTS are SMILES
    Format::SMILES,
    Format::EXTENDED_SMILES,
    Format::SMARTS,
    Format::EXTENDED_SMARTS,
    Format::HELM,
};

} // namespace

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
    ui->status_label->setText(text);
    ui->status_label->setStyleSheet(style);

    // load data into the format combo box
    auto interface_type = m_sketcher_model->getInterfaceType();
    auto cur_mol_type = m_sketcher_model->getMoleculeType();
    // first figure out whether we want atomistic formats, monomeric formats, or
    // both
    auto allowed_mol_type = FormatMoleculeType::BOTH;
    if (interface_type == InterfaceType::ATOMISTIC ||
        (!replace && cur_mol_type == MoleculeType::ATOMISTIC)) {
        allowed_mol_type = FormatMoleculeType::ATOMISTIC;
    } else if (interface_type == InterfaceType::MONOMERIC ||
               (!replace && cur_mol_type == MoleculeType::MONOMERIC)) {
        allowed_mol_type = FormatMoleculeType::MONOMERIC;
    }
    // then iterate through the available formats and load tha applicable ones
    // into the combo box
    ui->format_combo->clear();
    for (const auto& [format_name, format, format_mol_type] : FORMATS) {
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
        rdkit_extensions::to_rdkit(text, format);
        return true;
    } catch (const std::exception&) {
    }

    try {
        rdkit_extensions::to_rdkit_reaction(text, format);
        return true;
    } catch (const std::exception&) {
    }

    return false;
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
        for (const auto cur_format : AUTO_DETECT_FORMATS) {
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
    } else if (m_autodetected_format != std::nullopt) {
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
