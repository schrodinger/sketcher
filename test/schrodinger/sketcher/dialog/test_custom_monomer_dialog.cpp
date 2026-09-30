#define BOOST_TEST_MODULE Test_Sketcher

#include <string>
#include <utility>
#include <vector>

#include <QComboBox>
#include <QCoreApplication>
#include <QDialogButtonBox>
#include <QPushButton>
#include <QTextEdit>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

using rdkit_extensions::ChainType;

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
    CustomMonomerDialog dialog(ChainType::PEPTIDE);
    auto* sketcher = dialog.findChild<SketcherWidget*>();
    auto* ok_button = get_ok_button(dialog);
    BOOST_REQUIRE(sketcher != nullptr);

    BOOST_TEST(dialog.windowTitle() == "Sketch Custom Peptide Monomer");
    BOOST_TEST(!ok_button->isEnabled());

    dialog.addSMILES("CC");
    BOOST_TEST(ok_button->isEnabled());

    bool monomer_accepted = false;
    auto accepted_chain_type = ChainType::CHEM;
    QObject::connect(&dialog, &CustomMonomerDialog::customMonomerAccepted,
                     &dialog,
                     [&monomer_accepted, &accepted_chain_type](
                         const std::string&, const ChainType type) {
                         monomer_accepted = true;
                         accepted_chain_type = type;
                     });
    dialog.accept();
    BOOST_TEST(monomer_accepted);
    BOOST_TEST(static_cast<int>(accepted_chain_type) ==
               static_cast<int>(ChainType::PEPTIDE));
}

BOOST_AUTO_TEST_CASE(custom_monomer_dialog_titles_reflect_chain_type)
{
    CustomMonomerDialog peptide_dialog(ChainType::PEPTIDE);
    CustomMonomerDialog nucleic_acid_dialog(ChainType::RNA);
    CustomMonomerDialog chem_dialog(ChainType::CHEM);

    BOOST_TEST(peptide_dialog.windowTitle() == "Sketch Custom Peptide Monomer");
    BOOST_TEST(nucleic_acid_dialog.windowTitle() ==
               "Sketch Custom Nucleic Acid Monomer");
    BOOST_TEST(chem_dialog.windowTitle() == "Sketch Custom Chem Monomer");
}

/**
 * Numbered attachment points must be unique. Report every duplicated number in
 * numerical order and do not accept the custom monomer.
 */
BOOST_AUTO_TEST_CASE(custom_monomer_dialog_rejects_duplicate_attachment_points)
{
    const std::vector<std::pair<std::string, QString>> test_cases = {
        {"*C* |$_R1;;_R1$|", "Multiple R1 attachment points found. All "
                             "attachment points must be unique."},
        {"*C(*)(*)C(*)(*)* |$_R3;;_R1;_R2;;_R3;_R1;_R2$|",
         "Multiple R1, R2, and R3 attachment points found. All attachment "
         "points must be unique."}};

    for (const auto& [smiles, expected_error] : test_cases) {
        CustomMonomerDialog dialog(ChainType::PEPTIDE);
        bool accepted = false;
        QObject::connect(&dialog, &CustomMonomerDialog::customMonomerAccepted,
                         [&accepted]() { accepted = true; });
        dialog.setRequiredAttachmentPoints({2});
        dialog.addSMILES(smiles);

        dialog.accept();

        BOOST_TEST(!accepted);
        auto* message_box_dialog = dialog.findChild<MessageBoxDialog*>();
        BOOST_REQUIRE(message_box_dialog != nullptr);
        auto* error_text =
            message_box_dialog->findChild<QTextEdit*>("text_edit");
        BOOST_REQUIRE(error_text != nullptr);
        BOOST_TEST(error_text->toPlainText() == expected_error);
    }
}

/** Unique numbered attachment points continue to be accepted. */
BOOST_AUTO_TEST_CASE(custom_monomer_dialog_accepts_unique_attachment_points)
{
    CustomMonomerDialog dialog(ChainType::PEPTIDE);
    bool accepted = false;
    QObject::connect(&dialog, &CustomMonomerDialog::customMonomerAccepted,
                     [&accepted]() { accepted = true; });
    dialog.addSMILES("*C* |$_R1;;_R2$|");

    dialog.accept();

    BOOST_TEST(accepted);
    BOOST_TEST(dialog.findChild<MessageBoxDialog*>() == nullptr);
}

BOOST_AUTO_TEST_CASE(
    custom_monomer_dialog_warns_before_removing_bound_attachment_points)
{
    const std::vector<std::pair<std::vector<int>, QString>> test_cases = {
        {{3},
         "R3 has been removed from this monomer but is currently bound. "
         "Continuing will remove this connection."},
        {{1, 2, 3, 4},
         "R3 and R4 have been removed from this monomer but are currently "
         "bound. Continuing will remove these connections."}};

    for (const auto& [required_attachment_points, expected_warning] :
         test_cases) {
        CustomMonomerDialog dialog(ChainType::PEPTIDE);
        bool accepted = false;
        QObject::connect(&dialog, &CustomMonomerDialog::customMonomerAccepted,
                         [&accepted]() { accepted = true; });
        dialog.setRequiredAttachmentPoints(required_attachment_points);
        dialog.addSMILES("*C* |$_R1;;_R2$|");

        dialog.accept();

        BOOST_TEST(!accepted);
        auto* warning_dialog = dialog.findChild<MessageBoxDialog*>();
        BOOST_REQUIRE(warning_dialog != nullptr);
        auto* warning_text = warning_dialog->findChild<QTextEdit*>("text_edit");
        BOOST_REQUIRE(warning_text != nullptr);
        BOOST_TEST(warning_text->toPlainText() == expected_warning);

        auto* button_box =
            warning_dialog->findChild<QDialogButtonBox*>("button_box");
        BOOST_REQUIRE(button_box != nullptr);
        button_box->button(QDialogButtonBox::Cancel)->click();
        BOOST_TEST(!accepted);
        QCoreApplication::processEvents();
    }
}

BOOST_AUTO_TEST_CASE(
    custom_monomer_dialog_continues_after_attachment_point_warning)
{
    CustomMonomerDialog dialog(ChainType::PEPTIDE);
    bool accepted = false;
    QObject::connect(&dialog, &CustomMonomerDialog::customMonomerAccepted,
                     [&accepted]() { accepted = true; });
    dialog.setRequiredAttachmentPoints({3});
    dialog.addSMILES("*C* |$_R1;;_R2$|");

    dialog.accept();

    BOOST_TEST(!accepted);
    auto* warning_dialog = dialog.findChild<MessageBoxDialog*>();
    BOOST_REQUIRE(warning_dialog != nullptr);
    auto* button_box =
        warning_dialog->findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    button_box->button(QDialogButtonBox::Ok)->click();
    BOOST_TEST(accepted);
}

} // namespace sketcher
} // namespace schrodinger
