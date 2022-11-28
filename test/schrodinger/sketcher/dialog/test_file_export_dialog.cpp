#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_file_export_dialog

#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_file_export_dialog.h"
#include "../test_common.h"
#include "../test_sketcherScene.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

class TestFileExportDialog : public FileExportDialog
{
  public:
    TestFileExportDialog(SketcherModel* model) : FileExportDialog(model){};
    using FileExportDialog::getFileContent;
    using FileExportDialog::getValidExtensions;
    using FileExportDialog::m_ui;

    void setComboFormat(Format format)
    {
        auto value = QVariant::fromValue(format);
        auto index = m_ui->format_combo->findData(value);
        m_ui->format_combo->setCurrentIndex(index);
    };

    Format getComboFormat()
    {
        return m_ui->format_combo->currentData().value<Format>();
    };
};

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

BOOST_AUTO_TEST_CASE(test_FileExportDialog_standard)
{
    testSketcherScene scene;
    auto model = scene.getModel();

    // Confirm default format and default filename
    TestFileExportDialog dlg(model);
    BOOST_TEST(!model->hasReaction());
    BOOST_TEST(dlg.getComboFormat() == Format::MDL_MOLV3000);
    BOOST_TEST(dlg.m_ui->filename_le->text().toStdString() == "structure");
    // There are 6 available structure export formats
    BOOST_TEST(dlg.m_ui->format_combo->count() == 6);
    auto exts = dlg.getValidExtensions();
    // And MDL has 4 available extensions
    BOOST_TEST(exts.count() == 4);
    BOOST_TEST(exts.contains(".mol"));

    // Change format to PDB and confirm available extensions
    dlg.setComboFormat(Format::PDB);
    BOOST_TEST(dlg.getComboFormat() == Format::PDB);
    exts = dlg.getValidExtensions();
    BOOST_TEST(exts.count() == 2);
    BOOST_TEST(exts.contains(".pdb"));
}

BOOST_AUTO_TEST_CASE(test_FileExportDialog_reaction)
{
    testSketcherScene scene;
    auto model = scene.getModel();

    // Create a reaction
    QPointF pos(0.0, 0.0);
    scene.addReactionArrowAt(pos);

    // Confirm default format for reactions and default filename
    TestFileExportDialog dlg(model);
    BOOST_TEST(model->hasReaction());
    BOOST_TEST(dlg.getComboFormat() == Format::MDL_MOLV3000);
    BOOST_TEST(dlg.m_ui->filename_le->text().toStdString() == "structure");
    // There are only 2 available reaction formats
    BOOST_TEST(dlg.m_ui->format_combo->count() == 2);
    auto exts = dlg.getValidExtensions();
    // And MDL has only 1 available extension
    BOOST_TEST(exts.count() == 1);
    BOOST_TEST(exts.contains(".rxn"));

    // Change format to SMILES and confirm available extension
    dlg.setComboFormat(Format::SMILES);
    BOOST_TEST(dlg.getComboFormat() == Format::SMILES);
    exts = dlg.getValidExtensions();
    BOOST_TEST(exts.count() == 1);
    BOOST_TEST(exts.contains(".rsmi"));
}
