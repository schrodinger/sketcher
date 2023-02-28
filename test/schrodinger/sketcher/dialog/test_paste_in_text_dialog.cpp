#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>
#include "../test_common.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/dialog/paste_in_text_dialog.h"
#include "schrodinger/sketcher/ui/ui_paste_in_text_dialog.h"
#include <QSignalSpy>

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

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
        return ui->status_label->text().toStdString();
    }
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
    BOOST_TEST(args.at(0).toString().toStdString() == exp_text);
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

} // namespace sketcher
} // namespace schrodinger
