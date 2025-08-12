
#define BOOST_TEST_MODULE test_file_export_dialog

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
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
    BOOST_TEST(dlg.m_ui->format_combo->count() == 8); // export formats
    auto exts = dlg.getValidExtensions();
    BOOST_TEST(exts.size() == 8); // MDL available extensions
    BOOST_TEST(contains(exts, ".mol"));

    // Change format to PDB and confirm available extensions
    dlg.setComboFormat(Format::PDB);
    BOOST_TEST(dlg.getComboFormat() == Format::PDB);
    exts = dlg.getValidExtensions();
    BOOST_TEST(exts.size() == 6);
    BOOST_TEST(contains(exts, ".pdb"));
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
    // There are only 2 available reaction formats
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
