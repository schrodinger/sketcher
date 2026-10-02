#define BOOST_TEST_MODULE test_monomer_tool_widget
#include <boost/test/unit_test.hpp>

#include <QAbstractButton>
#include <QGridLayout>
#include <QPointer>
#include <QScrollArea>
#include <QScrollBar>
#include <QSignalSpy>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/amino_acid_symbol_popup.h"
#include "schrodinger/sketcher/widget/modular_tool_button.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"
#include "schrodinger/sketcher/widget/nucleic_acid_symbol_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Keep the existing column thresholds and make overflow buttons reachable by
 * scrolling in both kinds of monomer popup.
 */
BOOST_AUTO_TEST_CASE(monomer_popup_scrolling)
{
    for (const bool amino_acid : {true, false}) {
        for (const int count : {1, 20, 21, 80, 81, 200}) {
            BOOST_TEST_CONTEXT("amino_acid=" << amino_acid
                                             << ", count=" << count)
            {
                std::vector<rdkit_extensions::MonomerInfo> analogs(count - 1);
                for (int i = 0; i < count - 1; ++i) {
                    analogs[i].symbol = "M" + std::to_string(i);
                    analogs[i].name = "Monomer " + std::to_string(i);
                }
                std::unique_ptr<ModularPopup> popup;
                if (amino_acid) {
                    popup = std::make_unique<AminoAcidSymbolPopup>(
                        "A", "Alanine", analogs);
                } else {
                    popup = std::make_unique<NucleicAcidSymbolPopup>(
                        "A", "Adenine", analogs);
                }
                popup->show();
                QApplication::processEvents();
                auto* scroll_area = popup->findChild<QScrollArea*>();
                BOOST_REQUIRE((scroll_area != nullptr) == (count > 80));
                auto* grid = qobject_cast<QGridLayout*>(
                    scroll_area ? scroll_area->widget()->layout()
                                : popup->layout());
                BOOST_REQUIRE(grid != nullptr);
                const int columns = count <= 20 ? 4 : 8;
                BOOST_TEST(grid->columnCount() == std::min(count, columns));
                BOOST_TEST(grid->rowCount() == (count + columns - 1) / columns);
                const auto packets = popup->getButtonPackets();
                BOOST_REQUIRE_EQUAL(packets.size(), count);
                auto* last_button = packets.back().button;
                if (scroll_area != nullptr) {
                    auto* scrollbar = scroll_area->verticalScrollBar();
                    BOOST_TEST(scrollbar->isVisible());
                    BOOST_TEST(scrollbar->maximum() > 0);
                    BOOST_TEST(scroll_area->horizontalScrollBar()->maximum() ==
                               0);
                    BOOST_TEST(scroll_area->widget()->width() <=
                               scroll_area->viewport()->width());
                    scrollbar->setValue(scrollbar->maximum());
                    QApplication::processEvents();
                    const auto position = last_button->mapTo(
                        scroll_area->viewport(), QPoint(0, 0));
                    BOOST_TEST(position.y() >= 0);
                    BOOST_TEST(position.y() + last_button->height() <=
                               scroll_area->viewport()->height());
                }
                QSignalSpy spy(popup.get(), &ModularPopup::selectionChanged);
                last_button->click();
                BOOST_REQUIRE_EQUAL(spy.size(), 1);
                BOOST_TEST(spy.front().front().toInt() == count - 1);
                BOOST_TEST(!popup->isVisible());
            }
        }
    }
}

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
    auto second_widget = std::make_unique<MonomerToolWidget>();

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
    BOOST_REQUIRE(second_widget->findChild<QAbstractButton*>(
                      "analog_testAA_btn") != nullptr);
    // Widgets constructed after a load also see the current database.
    {
        MonomerToolWidget late_widget;
        BOOST_REQUIRE(late_widget.findChild<QAbstractButton*>(
                          "analog_testAA_btn") != nullptr);
    }
    // Subsequent updates must be safe after a subscriber is destroyed.
    second_widget.reset();
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
    BOOST_TEST(old_popup.isNull());
    BOOST_TEST(button->text() == "A");
    BOOST_TEST(widget.findChild<QAbstractButton*>("analog_testAA_btn") ==
               nullptr);
    BOOST_TEST(widget.findChild<QAbstractButton*>("na_analog_testNA_btn") ==
               nullptr);
    // Repeated refreshes must also replace existing core analog popups without
    // throwing an exception
    db.resetMonomerDefinitions();

    // Insertion into an existing database also refreshes the controls.
    db.loadMonomersFromJson("[]");
    db.insertMonomersFromJson(json);
    BOOST_REQUIRE(widget.findChild<QAbstractButton*>("analog_testAA_btn") !=
                  nullptr);
    QPointer<QWidget> inserted_popup = button->getPopupWidget();
    BOOST_CHECK_THROW(db.loadMonomersFromJson("invalid JSON"), std::exception);
    BOOST_TEST(button->getPopupWidget() == inserted_popup.data());
}

/** Unclassified buttons have no default and remember an explicit selection. */
BOOST_AUTO_TEST_CASE(unclassified_monomers)
{
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    db.loadMonomersFromJson("[]");
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    auto* aa_button =
        widget.findChild<ModularToolButton*>("aa_unclassified_btn");
    auto* na_button =
        widget.findChild<ModularToolButton*>("na_unclassified_btn");
    BOOST_REQUIRE(aa_button != nullptr);
    BOOST_REQUIRE(na_button != nullptr);
    const auto core_aa_count =
        db.getMonomersByNaturalAnalog(rdkit_extensions::ChainType::PEPTIDE)["X"]
            .size();
    BOOST_TEST(aa_button->isHidden() == (core_aa_count == 0));
    BOOST_TEST(na_button->isHidden());

    auto result = db.loadMonomersFromJson(R"([
        {"symbol":"testX","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CC","name":"Unclassified peptide","monomer_type":"backbone",
         "author":"test"},
        {"symbol":"testX2","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CN","name":"Another peptide","monomer_type":"backbone",
         "author":"test"},
        {"symbol":"testN","polymer_type":"RNA","natural_analog":"N",
         "smiles":"CCC","name":"Unclassified base","monomer_type":"branch",
         "author":"test"},
        {"symbol":"testN2","polymer_type":"RNA","natural_analog":"N",
         "smiles":"CCN","name":"Another base","monomer_type":"branch",
         "author":"test"}
    ])");
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), 4);
    for (bool amino_acid : {true, false}) {
        auto* button = amino_acid ? aa_button : na_button;
        const auto symbol = amino_acid ? "testX2" : "testN2";
        widget
            .findChild<QAbstractButton*>(amino_acid ? "amino_monomer_btn"
                                                    : "nucleic_monomer_btn")
            ->click();
        BOOST_TEST(!button->isHidden());
        BOOST_TEST(button->getEnumItem() == -1);
        BOOST_TEST(button->text() == "Unclassified");
        auto* popup = dynamic_cast<ModularPopup*>(button->getPopupWidget());
        BOOST_REQUIRE(popup != nullptr);
        const auto packets = popup->getButtonPackets();
        BOOST_REQUIRE_EQUAL(packets.size(),
                            (amino_acid ? core_aa_count : 0) + 2);
        auto* first_analog = popup->findChild<QToolButton*>(
            amino_acid ? "analog_testX_btn" : "na_analog_testN_btn");
        auto* second_analog = popup->findChild<QToolButton*>(
            amino_acid ? "analog_testX2_btn" : "na_analog_testN2_btn");
        BOOST_REQUIRE(first_analog != nullptr);
        BOOST_REQUIRE(second_analog != nullptr);
        BOOST_TEST(popup->findChild<QToolButton*>(
                       amino_acid ? "analog_X_btn" : "na_analog_N_btn") ==
                   nullptr);

        const auto key = amino_acid ? ModelKey::AMINO_ACID_TOOL
                                    : ModelKey::NUCLEIC_ACID_TOOL;
        const auto previous_tool = model.getValue(key);
        button->click();
        BOOST_TEST(popup->isVisible());
        BOOST_CHECK(model.getValue(key) == previous_tool);
        BOOST_TEST(!button->isChecked());
        popup->close();
        button->click();
        BOOST_TEST(popup->isVisible());
        first_analog->click();
        BOOST_TEST(!popup->isVisible());
        BOOST_TEST(button->isChecked());
        const auto first_id = button->getEnumItem();
        BOOST_TEST(first_id >= 0);
        BOOST_TEST(button->text() == "Unclassified");
        BOOST_TEST(button->styleSheet().contains(CUSTOM_MONOMER_BUTTON_STYLE));

        button->click();
        BOOST_TEST(popup->isVisible());
        second_analog->click();
        BOOST_TEST(!popup->isVisible());
        BOOST_TEST(button->getEnumItem() != first_id);
        BOOST_TEST(button->text() == "Unclassified");
        BOOST_TEST(button->styleSheet().contains(CUSTOM_MONOMER_BUTTON_STYLE));

        // Activate another button, then reuse the remembered analog.
        widget.findChild<QAbstractButton*>(amino_acid ? "ala_btn" : "na_a_btn")
            ->click();
        button->click();
        BOOST_TEST(!popup->isVisible());
        BOOST_TEST(button->isChecked());
        if (amino_acid) {
            BOOST_CHECK(model.getAminoAcidTool() ==
                        AminoAcidTool::UNCLASSIFIED);
            BOOST_TEST(model.getValueString(ModelKey::AMINO_ACID_SYMBOL) ==
                       symbol);
        } else {
            BOOST_CHECK(model.getNucleicAcidTool() ==
                        NucleicAcidTool::UNCLASSIFIED);
            const auto mutation = model.getValue(ModelKey::NUCLEIC_ACID_SYMBOL)
                                      .value<NucleicAcidMutation>();
            BOOST_CHECK(mutation.tool == NucleicAcidTool::UNCLASSIFIED);
            BOOST_TEST(mutation.symbol == symbol);
        }
    }

    QPointer<QWidget> old_aa_popup = aa_button->getPopupWidget();
    QPointer<QWidget> old_na_popup = na_button->getPopupWidget();
    db.loadMonomersFromJson("[]");
    BOOST_TEST(old_aa_popup.isNull());
    BOOST_TEST(old_na_popup.isNull());
    BOOST_TEST(aa_button->isHidden() == (core_aa_count == 0));
    BOOST_TEST(na_button->isHidden());
    BOOST_TEST(aa_button->getEnumItem() == -1);
    BOOST_TEST(na_button->getEnumItem() == -1);
    BOOST_TEST(aa_button->text() == "Unclassified");
    BOOST_TEST(na_button->text() == "Unclassified");
}

} // namespace sketcher
} // namespace schrodinger
