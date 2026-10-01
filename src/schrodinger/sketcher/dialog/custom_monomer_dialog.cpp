#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"

#include <algorithm>
#include <map>
#include <unordered_set>

#include <QDialogButtonBox>
#include <QPushButton>
#include <QStringList>

#include <stdexcept>

#include <rdkit/GraphMol/MolOps.h>

#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui/ui_custom_monomer_dialog.h"

using schrodinger::rdkit_extensions::ChainType;
using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::get_r_group_number;

namespace schrodinger
{
namespace sketcher
{

/**
 * @return a list of all R-groups that appear more than once in the specified
 * molecule
 */
static std::vector<int> get_duplicate_r_groups(const RDKit::ROMol& mol)
{
    std::map<unsigned int, unsigned int> r_group_counts;
    for (const auto* atom : mol.atoms()) {
        if (const auto r_group_num = get_r_group_number(atom)) {
            ++r_group_counts[*r_group_num];
        }
    }

    std::vector<int> duplicate_r_groups;
    for (const auto& [r_group_num, count] : r_group_counts) {
        if (count > 1) {
            duplicate_r_groups.push_back(static_cast<int>(r_group_num));
        }
    }
    return duplicate_r_groups;
}

/**
 * @return formatted text listing all R-groups in the given list of R-groups
 */
static QString format_r_group_list(const std::vector<int>& r_group_numbers)
{
    if (r_group_numbers.empty()) {
        return {};
    }

    QStringList r_groups;
    for (const auto r_group_num : r_group_numbers) {
        r_groups.append("R" + QString::number(r_group_num));
    }

    if (r_groups.size() == 1) {
        return r_groups.front();
    }

    const auto last_r_group = r_groups.takeLast();
    const auto separator = r_groups.size() == 1 ? " " : ", ";
    return r_groups.join(", ") + separator + "and " + last_r_group;
}

static QString chain_type_display_name(const ChainType chain_type)
{
    switch (chain_type) {
        case ChainType::PEPTIDE:
            return "Peptide";
        case ChainType::RNA:
            return "Nucleic Acid";
        case ChainType::CHEM:
            return "Chem";
        default:
            throw std::invalid_argument(
                "Custom monomer dialog does not support this chain type");
    }
}

CustomMonomerDialog::CustomMonomerDialog(const ChainType chain_type,
                                         QWidget* parent) :
    ModalDialog(parent),
    m_chain_type(chain_type)
{
    ui.reset(new Ui::CustomMonomerDialog());
    setupDialogUI(*ui);
    setStyleSheet(CUSTOM_MONOMER_DIALOG_STYLE);
    setWindowTitle("Sketch Custom " + chain_type_display_name(chain_type) +
                   " Monomer");
    ui->sketcher_widget->setInterfaceType(InterfaceType::ATOMISTIC);

    // remove the standard margins set by ModalDialog so that there's no gap
    // between the SketcherWidget and the edge of the dialog
    m_dlg_layout->setContentsMargins(0, 0, 0, 0);
    qobject_cast<QVBoxLayout*>(layout())->setContentsMargins(0, 0, 0, 0);

    connect(ui->sketcher_widget, &SketcherWidget::moleculeChanged, this,
            &CustomMonomerDialog::updateOkButton);
    connect(ui->sketcher_widget, &SketcherWidget::representationChanged, this,
            &CustomMonomerDialog::updateOkButton);
    updateOkButton();
}

CustomMonomerDialog::~CustomMonomerDialog() = default;

void CustomMonomerDialog::setRequiredAttachmentPoints(
    std::vector<int> required_attachment_points)
{
    std::ranges::sort(required_attachment_points);
    m_required_attachment_points.clear();
    std::ranges::unique_copy(required_attachment_points,
                             std::back_inserter(m_required_attachment_points));
}

void CustomMonomerDialog::addSMILES(const std::string& smiles)
{
    ui->sketcher_widget->addFromString(smiles, Format::EXTENDED_SMILES);
}

void CustomMonomerDialog::updateOkButton()
{
    bool valid = false;
    try {
        if (!ui->sketcher_widget->isEmpty()) {
            auto mol = ui->sketcher_widget->getRDKitMolecule();
            valid = mol->getNumAtoms() > 0 &&
                    RDKit::MolOps::getMolFrags(*mol, false).size() == 1;
            if (valid) {
                ui->sketcher_widget->getString(Format::EXTENDED_SMILES);
            }
        }
    } catch (const std::exception&) {
        valid = false;
    }
    ui->button_box->button(QDialogButtonBox::Ok)->setEnabled(valid);
}

/**
 * @return the test of the warning to use when the user has deleted the
 * specified required attachment points (i.e. attachment points that have a
 * bound connection)
 */
static QString
get_warning_text(const std::vector<int>& missing_attachment_points)
{
    const bool plural = missing_attachment_points.size() != 1;
    const auto attachment_points =
        format_r_group_list(missing_attachment_points);
    return attachment_points + (plural ? " have" : " has") +
           " been removed from this monomer but " + (plural ? "are" : "is") +
           " currently bound. Continuing will remove " +
           (plural ? "these connections." : "this connection.");
}

void CustomMonomerDialog::accept()
{
    const auto mol = ui->sketcher_widget->getRDKitMolecule();
    const auto duplicate_r_groups = get_duplicate_r_groups(*mol);
    if (!duplicate_r_groups.empty()) {
        show_error_dialog("Invalid Attachment Points",
                          "Multiple " +
                              format_r_group_list(duplicate_r_groups) +
                              " attachment points found. All attachment points "
                              "must be unique.",
                          this);
        return;
    }

    // warn before accepting if the user has deleted any required attachment
    // points (i.e. attachment points that have a bound connection)
    const auto missing_attachment_points =
        get_missing_required_attachment_points(*mol,
                                               m_required_attachment_points);
    const auto smiles = ui->sketcher_widget->getString(Format::EXTENDED_SMILES);
    if (!missing_attachment_points.empty()) {
        auto warning_text = get_warning_text(missing_attachment_points);
        auto* warning_dialog = show_warning_dialog("Remove Bound Connections?",
                                                   warning_text, this);
        connect(warning_dialog, &MessageBoxDialog::accepted, this,
                [this, smiles]() {
                    emit customMonomerAccepted(smiles, m_chain_type);
                    ModalDialog::accept();
                });
        return;
    }

    emit customMonomerAccepted(smiles, m_chain_type);
    ModalDialog::accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/custom_monomer_dialog.moc"
