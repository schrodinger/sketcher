#define BOOST_TEST_MODULE Test_Sketcher
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/paste_in_text_dialog.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_paste_in_text_dialog.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

using schrodinger::rdkit_extensions::Format;

namespace bdata = boost::unit_test::data;

namespace schrodinger
{
namespace sketcher
{

class TestPasteInTextDialog : public PasteInTextDialog
{
  public:
    TestPasteInTextDialog(SketcherModel* model) : PasteInTextDialog(model)
    {
    }
    std::string getText() const
    {
        return ui->structure_text_edit->toPlainText().toStdString();
    }

    void setText(const std::string& text)
    {
        ui->structure_text_edit->setPlainText(QString::fromStdString(text));
    }

    std::string getLabelText()
    {
        return ui->status_lbl->text().toStdString();
    }

    using PasteInTextDialog::ui;
};

/**
 * Verify that the dialog emits a signal with the text box text when accepted.
 */
BOOST_AUTO_TEST_CASE(text)
{

    SketcherModel model;
    TestPasteInTextDialog dlg(&model);
    QSignalSpy spy(&dlg, &TestPasteInTextDialog::textAccepted);
    std::string exp_text = "CCC";
    dlg.setText(exp_text);
    BOOST_TEST(dlg.getText() == exp_text);
    dlg.accept();

    BOOST_TEST(spy.count() == 1);
    auto args = spy.takeLast();
    BOOST_TEST(args.at(0).value<std::string>() == exp_text);
}

/**
 * Verify that the dialog label text matches the model state.
 */
BOOST_AUTO_TEST_CASE(label)
{

    SketcherModel model;
    TestPasteInTextDialog dlg(&model);

    // Verify initial state is as expected
    std::string replace_text =
        "Specified structure will <b>replace</b> Sketcher content";
    std::string add_text =
        "Specified structure will be added to Sketcher content";
    BOOST_TEST(model.getValueBool(ModelKey::NEW_STRUCTURES_REPLACE_CONTENT));
    BOOST_TEST(dlg.getLabelText() == replace_text);

    for (bool replace : {false, true}) {
        model.setValue(ModelKey::NEW_STRUCTURES_REPLACE_CONTENT, replace);
        auto exp_text = replace ? replace_text : add_text;
        BOOST_TEST(dlg.getLabelText() == exp_text);
    }
}

/**
 * Make sure that the format combo box gets loaded with the correct number of
 * formats for atomistic vs monomeric vs both. Also make sure that autodetect
 * only detects acceptable formats.
 */
BOOST_DATA_TEST_CASE(
    test_format_combo_contents,
    bdata::make(std::vector<std::tuple<InterfaceTypeType, MoleculeType, bool,
                                       int, bool, bool>>{
        {InterfaceType::ATOMISTIC, MoleculeType::EMPTY, false, 13, true, false},
        {InterfaceType::ATOMISTIC, MoleculeType::ATOMISTIC, false, 13, true,
         false},
        {InterfaceType::ATOMISTIC, MoleculeType::EMPTY, true, 13, true, false},
        {InterfaceType::ATOMISTIC, MoleculeType::ATOMISTIC, true, 13, true,
         false},
        {InterfaceType::MONOMERIC, MoleculeType::EMPTY, false, 5, false, true},
        {InterfaceType::MONOMERIC, MoleculeType::MONOMERIC, false, 5, false,
         true},
        {InterfaceType::MONOMERIC, MoleculeType::EMPTY, true, 5, false, true},
        {InterfaceType::MONOMERIC, MoleculeType::MONOMERIC, true, 5, false,
         true},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::EMPTY, false, 17,
         true, true},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::ATOMISTIC, false,
         13, true, false},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::MONOMERIC, false,
         5, false, true},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::EMPTY, true, 17,
         true, true},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::ATOMISTIC, true,
         17, true, true},
        {InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::MONOMERIC, true,
         17, true, true},
    }),
    interface_type, molecule_type, replace_contents, num_formats,
    accepts_atomistic, accepts_monomeric)
{
    SketcherModel model;
    model.setValues(
        {{ModelKey::INTERFACE_TYPE, interface_type},
         {ModelKey::MOLECULE_TYPE, QVariant::fromValue(molecule_type)},
         {ModelKey::NEW_STRUCTURES_REPLACE_CONTENT, replace_contents}});
    TestPasteInTextDialog dlg(&model);
    BOOST_TEST(dlg.ui->format_combo->count() == num_formats);
    BOOST_TEST(dlg.ui->format_combo->currentText() == "Autodetect");

    dlg.setText("CCC");
    QString exp_text =
        accepts_atomistic ? "Autodetect (SMILES)" : "Autodetect (unrecognized)";
    BOOST_TEST(dlg.ui->format_combo->currentText() == exp_text);

    dlg.setText("[#6]");
    exp_text =
        accepts_atomistic ? "Autodetect (SMARTS)" : "Autodetect (unrecognized)";
    BOOST_TEST(dlg.ui->format_combo->currentText() == exp_text);

    dlg.setText("PEPTIDE1{K.W.L.N.A.L.L.H.H.G.L}$$$$V2.0");
    exp_text =
        accepts_monomeric ? "Autodetect (HELM)" : "Autodetect (unrecognized)";
    BOOST_TEST(dlg.ui->format_combo->currentText() == exp_text);

    // gibberish text should always be reported as unrecognized
    dlg.setText("](@*&|;");
    BOOST_TEST(dlg.ui->format_combo->currentText() ==
               "Autodetect (unrecognized)");

    // clearing the text should remove the parenthetical
    dlg.setText("");
    BOOST_TEST(dlg.ui->format_combo->currentText() == "Autodetect");
}

/**
 * Make sure that the list of formats includes all of the reaction formats.
 * These aren't explicitly loaded into the dialog, as we assume that they're a
 * subset of the atomistic formats. This test ensures that this assumption is
 * still valid.
 */
BOOST_AUTO_TEST_CASE(test_reaction_formats)
{
    SketcherModel model;
    model.setValues(
        {{ModelKey::INTERFACE_TYPE, InterfaceType::ATOMISTIC_OR_MONOMERIC}});
    TestPasteInTextDialog dlg(&model);
    for (auto& [rxn_format, format_name] : get_reaction_import_formats()) {
        auto found_idx =
            dlg.ui->format_combo->findData(QVariant::fromValue(rxn_format));
        BOOST_TEST(found_idx > 0);
    }
}

} // namespace sketcher
} // namespace schrodinger
