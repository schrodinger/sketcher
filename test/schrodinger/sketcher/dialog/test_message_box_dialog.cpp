#define BOOST_TEST_MODULE Test_Sketcher

#include <QDialogButtonBox>
#include <QWidget>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_show_warning_dialog)
{
    QWidget parent;

    show_warning_dialog("Warning", "Warning text", &parent);

    auto* dialog = parent.findChild<MessageBoxDialog*>();
    BOOST_REQUIRE(dialog != nullptr);
    auto* button_box = dialog->findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    BOOST_TEST(button_box->standardButtons().testFlag(QDialogButtonBox::Ok));
    BOOST_TEST(
        button_box->standardButtons().testFlag(QDialogButtonBox::Cancel));
}

BOOST_AUTO_TEST_CASE(test_show_error_dialog)
{
    QWidget parent;

    show_error_dialog("Error", "Error text", &parent);

    auto* dialog = parent.findChild<MessageBoxDialog*>();
    BOOST_REQUIRE(dialog != nullptr);
    auto* button_box = dialog->findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    BOOST_TEST(button_box->standardButtons().testFlag(QDialogButtonBox::Ok));
    BOOST_TEST(
        !button_box->standardButtons().testFlag(QDialogButtonBox::Cancel));
}

BOOST_AUTO_TEST_CASE(test_show_information_dialog)
{
    QWidget parent;

    show_information_dialog("Information", "Information text", &parent);

    auto* dialog = parent.findChild<MessageBoxDialog*>();
    BOOST_REQUIRE(dialog != nullptr);
    auto* button_box = dialog->findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    BOOST_TEST(button_box->standardButtons().testFlag(QDialogButtonBox::Ok));
    BOOST_TEST(
        !button_box->standardButtons().testFlag(QDialogButtonBox::Cancel));
}

} // namespace sketcher
} // namespace schrodinger
