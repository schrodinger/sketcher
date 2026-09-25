#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"

#include <QDialogButtonBox>
#include <QPushButton>

#include <stdexcept>

#include <rdkit/GraphMol/MolOps.h>

#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui/ui_custom_monomer_dialog.h"

using schrodinger::rdkit_extensions::ChainType;
using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

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

void CustomMonomerDialog::accept()
{
    auto smiles = ui->sketcher_widget->getString(Format::EXTENDED_SMILES);
    emit customMonomerAccepted(smiles, m_chain_type);
    ModalDialog::accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/custom_monomer_dialog.moc"
