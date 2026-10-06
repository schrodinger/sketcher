#define BOOST_TEST_MODULE Test_Sketcher

#include <string>
#include <utility>
#include <vector>

#include <QComboBox>
#include <QCoreApplication>
#include <QDialogButtonBox>
#include <QLayout>
#include <QPushButton>
#include <QTextEdit>
#include <QToolButton>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/widget/sketcher_side_bar.h"
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
 * Resizing switches footer placement without raising the compact minimum or
 * changing the threshold. The footer remains visible and its signals work
 * after repeated reparenting.
 */
static void check_footer_placement(CustomMonomerDialog& dialog)
{
    dialog.setAttribute(Qt::WA_DeleteOnClose, false);
    auto* sketcher = dialog.findChild<SketcherWidget*>();
    auto* footer = dialog.findChild<QWidget*>("button_bar");
    auto* view = dialog.findChild<QWidget*>("view");
    auto* sidebar = dialog.findChild<QWidget*>("side_bar_wdg");
    auto* button_box = dialog.findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(sketcher != nullptr);
    BOOST_REQUIRE(footer != nullptr);
    BOOST_REQUIRE(view != nullptr);
    BOOST_REQUIRE(sidebar != nullptr);
    BOOST_REQUIRE(button_box != nullptr);

    dialog.resize(dialog.sizeHint() + QSize(100, 100));
    dialog.show();
    QCoreApplication::processEvents();
    const auto compact_minimum = dialog.minimumSizeHint();
    const int full_footer_height = dialog.layout()->totalMinimumSize().height();
    BOOST_TEST(full_footer_height > compact_minimum.height());
    BOOST_TEST(footer->parentWidget() != sketcher);
    BOOST_TEST(footer->mapTo(&dialog, QPoint()).y() >=
               sidebar->mapTo(&dialog, QPoint()).y() + sidebar->height());

    for (int i = 0; i < 3; ++i) {
        dialog.resize(dialog.width(), full_footer_height - 1);
        QCoreApplication::processEvents();
        BOOST_TEST(dialog.height() == full_footer_height - 1);
        BOOST_TEST(footer->parentWidget() == sketcher);
        BOOST_TEST(footer->isVisible());
        BOOST_TEST(footer->mapTo(&dialog, QPoint()).y() >=
                   view->mapTo(&dialog, QPoint()).y() + view->height());
        BOOST_TEST(dialog.minimumSizeHint().height() ==
                   compact_minimum.height());
        BOOST_TEST(dialog.minimumSizeHint().width() == compact_minimum.width());
        BOOST_TEST(dialog.minimumHeight() == compact_minimum.height());
        BOOST_TEST(dialog.layout()->totalMinimumSize().height() ==
                   compact_minimum.height());

        dialog.resize(dialog.width(), full_footer_height);
        QCoreApplication::processEvents();
        BOOST_TEST(footer->parentWidget() != sketcher);
        BOOST_TEST(footer->isVisible());
        BOOST_TEST(dialog.minimumSizeHint().height() ==
                   compact_minimum.height());
        BOOST_TEST(dialog.layout()->totalMinimumSize().height() ==
                   full_footer_height);
    }

    dialog.resize(compact_minimum);
    QCoreApplication::processEvents();
    BOOST_TEST(dialog.height() == compact_minimum.height());
    BOOST_TEST(dialog.width() == compact_minimum.width());
    BOOST_TEST(footer->parentWidget() == sketcher);
    BOOST_TEST(view->height() >= view->minimumSizeHint().height());
    BOOST_TEST(button_box->width() >= button_box->minimumSizeHint().width());

    bool rejected = false;
    QObject::connect(&dialog, &QDialog::rejected,
                     [&rejected]() { rejected = true; });
    button_box->button(QDialogButtonBox::Cancel)->click();
    BOOST_TEST(rejected);
}

BOOST_AUTO_TEST_CASE(custom_monomer_dialog_footer_follows_available_height)
{
    CustomMonomerDialog dialog(ChainType::PEPTIDE);
    check_footer_placement(dialog);
}

/**
 * Exercise title-bar and border accounting on native builds as well as WASM,
 * and make the View column determine the compact minimum height.
 */
BOOST_AUTO_TEST_CASE(
    custom_monomer_dialog_footer_with_title_bar_and_tall_footer)
{
    class DialogWithTitleBar : public CustomMonomerDialog
    {
      public:
        DialogWithTitleBar() : CustomMonomerDialog(ChainType::PEPTIDE)
        {
            if (m_title_bar == nullptr) {
                m_title_bar = new CustomTitleBar(windowTitle(), this);
                qobject_cast<QVBoxLayout*>(layout())->insertWidget(0,
                                                                   m_title_bar);
            }
            configureWasmTitleBar();
            setStyleSheet(
                styleSheet() +
                "QDialog { border: 1px solid #b5b5b5; }"
                "QDialogButtonBox > QPushButton { padding: 4px 20px; }");
            auto* sidebar = findChild<QWidget*>("side_bar_wdg");
            auto* footer = findChild<QWidget*>("button_bar");
            footer->setMinimumHeight(sidebar->minimumSizeHint().height());
        }
    } dialog;
    check_footer_placement(dialog);
}

/**
 * A title bar keeps its normal height when space permits, shrinks to the
 * polished toggle height minus the border margins at the compact minimum,
 * and grows again on resize.
 */
BOOST_AUTO_TEST_CASE(custom_monomer_dialog_title_bar_can_shrink)
{
    class DialogWithTitleBar : public CustomMonomerDialog
    {
      public:
        DialogWithTitleBar() : CustomMonomerDialog(ChainType::PEPTIDE)
        {
            if (m_title_bar == nullptr) {
                m_title_bar = new CustomTitleBar(windowTitle(), this);
                qobject_cast<QVBoxLayout*>(layout())->insertWidget(0,
                                                                   m_title_bar);
            }
            // Use a smaller toggle icon size so this checks an actual range
            // regardless of the platform's default widget style.
            findChild<QToolButton*>("atomistic_btn")
                ->setIconSize(QSize(16, 16));
            findChild<QToolButton*>("monomeric_btn")
                ->setIconSize(QSize(16, 16));
            configureWasmTitleBar();
        }

        CustomTitleBar* titleBar() const
        {
            return m_title_bar;
        }
    } dialog;
    dialog.setAttribute(Qt::WA_DeleteOnClose, false);
    auto* title_bar = dialog.titleBar();
    const int normal_height = title_bar->maximumHeight();
    const int minimum_height =
        dialog.findChild<SketcherSideBar*>()->getInterfaceToggleHeight() - 2;
    BOOST_REQUIRE(minimum_height < normal_height);
    BOOST_TEST(title_bar->minimumHeight() == minimum_height);

    dialog.resize(dialog.sizeHint() + QSize(100, 100));
    dialog.show();
    QCoreApplication::processEvents();
    BOOST_TEST(title_bar->height() == normal_height);
    const auto compact_minimum = dialog.minimumSizeHint();

    dialog.resize(compact_minimum);
    QCoreApplication::processEvents();
    BOOST_TEST(dialog.height() == compact_minimum.height());
    BOOST_TEST(title_bar->height() == minimum_height);
    BOOST_TEST(dialog.layout()->totalMinimumSize().height() ==
               compact_minimum.height());

    dialog.resize(dialog.sizeHint() + QSize(100, 100));
    QCoreApplication::processEvents();
    BOOST_TEST(title_bar->height() == normal_height);
    BOOST_TEST(dialog.minimumSizeHint().height() == compact_minimum.height());
    dialog.close();
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
