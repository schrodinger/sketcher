
#define BOOST_TEST_MODULE test_file_export_dialog

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_file_export_dialog.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

class TestFileExportDialog : public FileExportDialog
{
  public:
    TestFileExportDialog(SketcherModel* model) : FileExportDialog(model){};
    using FileExportDialog::getFileContent;
    using FileExportDialog::getValidExtensions;
    using FileExportDialog::m_current_format_list;
    using FileExportDialog::m_ui;

    void setComboFormat(Format format)
    {
        // Find the index in the format list that matches this format
        for (int i = 0; i < static_cast<int>(m_current_format_list.size());
             i++) {
            if (std::get<0>(m_current_format_list[i]) == format) {
                // Find the combo box item with this index
                for (int j = 0; j < m_ui->format_combo->count(); j++) {
                    if (m_ui->format_combo->itemData(j).toInt() == i) {
                        m_ui->format_combo->setCurrentIndex(j);
                        return;
                    }
                }
            }
        }
    };

    Format getComboFormat()
    {
        int format_index = m_ui->format_combo->currentData().toInt();
        return std::get<0>(m_current_format_list[format_index]);
    };
};

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

bool contains(const std::vector<std::string>& vec, const std::string& str)
{
    return std::find(vec.begin(), vec.end(), str) != vec.end();
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_standard)
{
    SketcherModel model;
    BOOST_TEST(!model.hasReaction());

    // Confirm default format and default filename
    TestFileExportDialog dlg(&model);
    BOOST_TEST(dlg.getComboFormat() == Format::MDL_MOLV3000);
    BOOST_TEST(dlg.m_ui->filename_le->text().toStdString() == "structure");
    // 12 atomistic entries (compressed variants listed separately) + HELM +
    // FASTA = 14
    BOOST_TEST(dlg.m_ui->format_combo->count() == 14);
    auto exts = dlg.getValidExtensions();
    BOOST_TEST(exts.size() == 4); // MDL uncompressed extensions only
    BOOST_TEST(contains(exts, ".mol"));

    // Change format to PDB and confirm available extensions
    dlg.setComboFormat(Format::PDB);
    BOOST_TEST(dlg.getComboFormat() == Format::PDB);
    exts = dlg.getValidExtensions();
    BOOST_TEST(exts.size() == 2); // PDB uncompressed extensions only
    BOOST_TEST(contains(exts, ".pdb"));

    // HELM and FASTA are always available regardless of model state;
    // selecting them produces sequence extensions.
    dlg.setComboFormat(Format::HELM);
    BOOST_TEST(contains(dlg.getValidExtensions(), ".helm"));
    dlg.setComboFormat(Format::FASTA);
    BOOST_TEST(contains(dlg.getValidExtensions(), ".fasta"));
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_reaction)
{
    SketcherModel model;
    QObject::connect(&model, &SketcherModel::reactionCountRequested,
                     []() { return 1; });
    BOOST_TEST(model.hasReaction());

    // Confirm default format for reactions and default filename
    TestFileExportDialog dlg(&model);

    BOOST_TEST(dlg.getComboFormat() == Format::MDL_MOLV3000);
    BOOST_TEST(dlg.m_ui->filename_le->text().toStdString() == "structure");
    // There are 2 reaction formats
    BOOST_TEST(dlg.m_ui->format_combo->count() == 2);
    auto exts = dlg.getValidExtensions();
    // And MDL has only 1 available extension
    BOOST_TEST(exts.size() == 1);
    BOOST_TEST(contains(exts, ".rxn"));

    // Change format to SMILES and confirm available extension
    dlg.setComboFormat(Format::SMILES);
    BOOST_TEST(dlg.getComboFormat() == Format::SMILES);
    exts = dlg.getValidExtensions();
    BOOST_TEST(exts.size() == 1);
    BOOST_TEST(contains(exts, ".rsmi"));
}

BOOST_AUTO_TEST_CASE(export_button)
{
    SketcherModel model;
    TestFileExportDialog dlg(&model);
#ifdef __EMSCRIPTEN__
    BOOST_TEST(dlg.m_ui->export_btn->text().toStdString() == "Download");
#else
    BOOST_TEST(dlg.m_ui->export_btn->text().toStdString() == "Save...");
#endif
}

// Wire the dialog to a real mol so its export plumbing (combo selection ->
// getFileContent() -> exportTextRequested signal) can be exercised end-to-end
// against rdkit_extensions::to_string.
static std::unique_ptr<TestFileExportDialog>
make_dialog_for_mol(SketcherModel& model,
                    const boost::shared_ptr<RDKit::RWMol>& mol)
{
    auto dlg = std::make_unique<TestFileExportDialog>(&model);
    QObject::connect(
        dlg.get(), &FileExportDialog::exportTextRequested,
        [mol](Format format) {
            return QString::fromStdString(
                schrodinger::rdkit_extensions::to_string(*mol, format));
        });
    return dlg;
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_helm_roundtrip)
{
    using schrodinger::rdkit_extensions::to_rdkit;
    auto source_mol = to_rdkit("PEPTIDE1{A.G.L}$$$$V2.0", Format::HELM);

    SketcherModel model;
    auto dlg = make_dialog_for_mol(model, source_mol);
    dlg->setComboFormat(Format::HELM);
    auto exported = dlg->getFileContent().toStdString();

    auto roundtripped = to_rdkit(exported, Format::HELM);
    BOOST_TEST(roundtripped->getNumAtoms() == source_mol->getNumAtoms());
    BOOST_TEST(roundtripped->getNumBonds() == source_mol->getNumBonds());
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_fasta_roundtrip)
{
    using schrodinger::rdkit_extensions::to_rdkit;
    auto source_mol = to_rdkit("PEPTIDE1{A.G.L}$$$$V2.0", Format::HELM);

    SketcherModel model;
    auto dlg = make_dialog_for_mol(model, source_mol);
    dlg->setComboFormat(Format::FASTA);
    auto exported = dlg->getFileContent().toStdString();

    // FASTA reads require the specific sub-format (peptide vs DNA vs RNA);
    // the exporter always emits Format::FASTA.
    auto roundtripped = to_rdkit(exported, Format::FASTA_PEPTIDE);
    BOOST_TEST(roundtripped->getNumAtoms() == source_mol->getNumAtoms());
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_monomeric_to_atomistic_export)
{
    using schrodinger::rdkit_extensions::to_rdkit;
    auto source_mol = to_rdkit("PEPTIDE1{A.G.L}$$$$V2.0", Format::HELM);
    BOOST_TEST(schrodinger::rdkit_extensions::isMonomeric(*source_mol));

    SketcherModel model;
    auto dlg = make_dialog_for_mol(model, source_mol);
    // Picking an atomistic format on a monomeric scene should expand the
    // monomers into their constituent atoms before serializing.
    dlg->setComboFormat(Format::MDL_MOLV3000);
    auto exported = dlg->getFileContent().toStdString();

    auto roundtripped = to_rdkit(exported, Format::MDL_MOLV3000);
    BOOST_TEST(!schrodinger::rdkit_extensions::isMonomeric(*roundtripped));
    BOOST_TEST(roundtripped->getNumAtoms() > source_mol->getNumAtoms());
}
