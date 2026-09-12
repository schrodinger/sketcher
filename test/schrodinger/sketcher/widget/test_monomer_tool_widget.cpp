#define BOOST_TEST_MODULE test_monomer_tool_widget
#include <boost/test/unit_test.hpp>

#include <QAbstractButton>
#include <QPointer>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/widget/modular_tool_button.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Verify that unknown monomer buttons use the unknown monomer styling. If these
 * buttons are converted to ModularToolButtons, then the ModularToolButton
 * styling will inadvertently overwrite the UNKNOWN_MONOMER_STYLE.
 */
BOOST_AUTO_TEST_CASE(unknown_monomer_button_styles)
{
    MonomerToolWidget widget;

    for (const auto* button_name : {"unk_btn", "na_n_btn"}) {
        const auto* button = widget.findChild<QAbstractButton*>(button_name);
        BOOST_REQUIRE(button != nullptr);
        const auto style_sheet = button->styleSheet();
        BOOST_TEST_CONTEXT(button->objectName().toStdString())
        {
            BOOST_TEST(style_sheet.contains(QStringLiteral("italic")));
            BOOST_TEST(style_sheet.contains(QStringLiteral("#606060")));
        }
    }
}

/**
 * Update the monomer database and confirm that the monomer button popups get
 * updated with the new monomer definitions
 */
BOOST_AUTO_TEST_CASE(refresh_monomer_popups)
{
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    // make sure that we undo any changes to the monomer database when the test
    // finishes
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    db.resetMonomerDefinitions();
    auto scene = TestScene::getScene();
    MonomerToolWidget widget;
    widget.setModel(scene->m_sketcher_model);

    // read in new monomers and confirm that they're added to the relevant pop
    // ups
    const auto json = R"([
        {"symbol":"testAA","polymer_type":"PEPTIDE","natural_analog":"A",
         "smiles":"CC","name":"Test amino acid","monomer_type":"backbone",
         "author":"test","pdbcode":"TAA"},
        {"symbol":"testNA","polymer_type":"RNA","natural_analog":"A",
         "smiles":"CCC","name":"Test nucleic acid","monomer_type":"branch",
         "author":"test","pdbcode":"TNA"}
    ])";
    auto result = db.loadMonomersFromJson(json);
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), 2);
    widget.updateMonomerButtons();
    for (const auto* name : {"ala_btn", "na_a_btn"}) {
        auto* button = widget.findChild<ModularToolButton*>(name);
        BOOST_REQUIRE(button != nullptr);
        auto* popup = button->getPopupWidget();
        BOOST_REQUIRE(popup != nullptr);
        auto* view = dynamic_cast<SketcherView*>(popup);
        BOOST_REQUIRE(view != nullptr);
        BOOST_TEST(view->getModel() == scene->m_sketcher_model);
        auto* analog = popup->findChild<QAbstractButton*>(
            QString(name) == "ala_btn" ? "analog_testAA_btn"
                                       : "na_analog_testNA_btn");
        BOOST_REQUIRE_MESSAGE(analog != nullptr, name);
        button->setEnumItem(1);
    }

    // reset the monomer definitions and confirm that the new monomers have been
    // removed from the pop ups
    auto* button = widget.findChild<ModularToolButton*>("ala_btn");
    QPointer<QWidget> old_popup = button->getPopupWidget();
    db.resetMonomerDefinitions();
    widget.updateMonomerButtons();
    BOOST_TEST(old_popup.isNull());
    BOOST_TEST(button->text() == "A");
    BOOST_TEST(widget.findChild<QAbstractButton*>("analog_testAA_btn") ==
               nullptr);
    BOOST_TEST(widget.findChild<QAbstractButton*>("na_analog_testNA_btn") ==
               nullptr);
    // Repeated refreshes must also replace existing core analog popups without
    // throwing an exception
    widget.updateMonomerButtons();
}

} // namespace sketcher
} // namespace schrodinger
