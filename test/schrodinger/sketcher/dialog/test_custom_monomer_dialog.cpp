#define BOOST_TEST_MODULE Test_Sketcher

#include <QDialogButtonBox>
#include <QPushButton>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"
#include "schrodinger/sketcher/sketcher_widget.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

static QPushButton* get_ok_button(const QWidget& dialog)
{
    auto* button_box = dialog.findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    return button_box->button(QDialogButtonBox::Ok);
}

/**
 * Make sure that the OK button is enabled only when the dialog is non-empty
 */
BOOST_AUTO_TEST_CASE(custom_monomer_dialog_validation_and_acceptance)
{
    CustomMonomerDialog dialog;
    auto* sketcher = dialog.findChild<SketcherWidget*>();
    auto* ok_button = get_ok_button(dialog);
    BOOST_REQUIRE(sketcher != nullptr);

    BOOST_TEST(dialog.windowTitle() == "Sketch Custom Monomer");
    BOOST_TEST(!ok_button->isEnabled());

    dialog.addSMILES("CC");
    BOOST_TEST(ok_button->isEnabled());
}

} // namespace sketcher
} // namespace schrodinger
