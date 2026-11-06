#define BOOST_TEST_MODULE test_monomer_tool_widget
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>

#include <QAbstractButton>
#include <QAction>
#include <QButtonGroup>
#include <QEnterEvent>
#include <QGridLayout>
#include <QGraphicsSceneMouseEvent>
#include <QMenu>
#include <QPointer>
#include <QProxyStyle>
#include <QScrollArea>
#include <QScrollBar>
#include <QSignalSpy>
#include <QStackedWidget>
#include <QStyleOptionToolButton>
#include <QTest>
#include <QTimer>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"
#include "schrodinger/sketcher/widget/amino_acid_symbol_popup.h"
#include "schrodinger/sketcher/widget/custom_monomer_button.h"
#include "schrodinger/sketcher/widget/modular_tool_button.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"
#include "schrodinger/sketcher/widget/nucleic_acid_symbol_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace rdkit_extensions
{
MAKE_ENUM_LOGGABLE(ChainType);
}

namespace sketcher
{

using rdkit_extensions::ChainType;
namespace bdata = boost::unit_test::data;

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
 * Clicking opens the custom monomer menu asynchronously without opening a
 * drawing dialog or changing the button's check state. Holding has no special
 * behavior.
 */
BOOST_AUTO_TEST_CASE(custom_monomer_menu_click)
{
    auto scene = TestScene::getScene();
    for (bool checked : {false, true}) {
        MonomerToolWidget widget;
        widget.setModel(scene->m_sketcher_model);
        widget.show();
        auto* button =
            widget.findChild<CustomMonomerButton*>("custom_monomer_btn");
        BOOST_REQUIRE(button != nullptr);
        auto* menu = button->findChild<QMenu*>();
        BOOST_REQUIRE(menu != nullptr);
        // Installing this as a built-in menu would reintroduce QMenu::exec().
        BOOST_TEST(button->menu() == nullptr);
        button->setChecked(checked);
        QSignalSpy clicked(button, &QToolButton::clicked);
        const auto draw_tool = scene->m_sketcher_model->getDrawTool();

        QTest::mousePress(button, Qt::LeftButton);
        QTest::qWait(300);
        BOOST_TEST(!menu->isVisible());
        BOOST_TEST(button->isChecked() == checked);
        BOOST_TEST(clicked.count() == 0);
        BOOST_TEST(widget.findChild<CustomMonomerDialog*>() == nullptr);

        // Ensure a regression to a blocking menu cannot hang this test.
        bool watchdog_fired = false;
        QTimer watchdog;
        watchdog.setSingleShot(true);
        QObject::connect(&watchdog, &QTimer::timeout, menu,
                         [menu, &watchdog_fired]() {
                             watchdog_fired = true;
                             menu->hide();
                         });
        watchdog.start(500);
        QTest::mouseRelease(button, Qt::LeftButton);
        watchdog.stop();
        BOOST_TEST(!watchdog_fired);
        BOOST_TEST(menu->isVisible());
        BOOST_TEST(clicked.count() == 1);
        BOOST_TEST(button->isChecked() == checked);
        BOOST_TEST(widget.findChild<CustomMonomerDialog*>() == nullptr);
        menu->hide();

        QEnterEvent enter(button->rect().center(), button->rect().center(),
                          button->mapToGlobal(button->rect().center()));
        QApplication::sendEvent(button, &enter);
        QTest::mouseClick(button, Qt::LeftButton);
        BOOST_TEST(menu->isVisible());
        BOOST_TEST(clicked.count() == 2);
        BOOST_TEST(button->isChecked() == checked);
        BOOST_TEST(widget.findChild<CustomMonomerDialog*>() == nullptr);
        BOOST_TEST(scene->m_sketcher_model->getDrawTool() == draw_tool);
        // Dismissal through the title bar or another window need not deliver
        // a leave event to the button that opened the menu.
        menu->hide();
        BOOST_TEST(!button->underMouse());
        BOOST_TEST(!button->isDown());
        BOOST_TEST(button->isChecked() == checked);

        menu->actions().front()->trigger();
        auto* dialog = button->findChild<CustomMonomerDialog*>();
        BOOST_REQUIRE(dialog != nullptr);
        BOOST_TEST(button->isChecked() == checked);
        dialog->reject();
        BOOST_TEST(button->isChecked() == checked);
        BOOST_TEST(scene->m_sketcher_model->getDrawTool() == draw_tool);
    }
}

/**
 * Menu actions follow the selected tab and forward accepted custom monomers
 * through the button to the model.
 */
BOOST_DATA_TEST_CASE(custom_monomer_menu_actions,
                     bdata::make(std::vector<ChainType>{
                         ChainType::PEPTIDE, ChainType::RNA, ChainType::CHEM}),
                     chain_type)
{
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    auto* button = widget.findChild<CustomMonomerButton*>("custom_monomer_btn");
    BOOST_REQUIRE(button != nullptr);
    auto* menu = button->findChild<QMenu*>();
    BOOST_REQUIRE(menu != nullptr);
    BOOST_REQUIRE_EQUAL(menu->actions().size(), 3);
    BOOST_TEST(menu->actions()[0]->text() == "Sketch Residue...");
    if (chain_type == ChainType::RNA) {
        widget.findChild<QToolButton*>("nucleic_monomer_btn")->click();
        BOOST_TEST(menu->actions()[0]->text() == "Sketch Base...");
        widget.findChild<QToolButton*>("amino_monomer_btn")->click();
        BOOST_TEST(menu->actions()[0]->text() == "Sketch Residue...");
        widget.findChild<QToolButton*>("nucleic_monomer_btn")->click();
    }
    menu->actions()[chain_type == ChainType::CHEM ? 2 : 0]->trigger();
    auto* dialog = button->findChild<CustomMonomerDialog*>();
    BOOST_REQUIRE(dialog != nullptr);
    BOOST_TEST(dialog->isVisible());
    BOOST_TEST(!button->isChecked());
    const QString expected_title =
        chain_type == ChainType::PEPTIDE ? "Define Custom Peptide Residue"
        : chain_type == ChainType::RNA   ? "Define Custom Nucleic Acid Base"
                                         : "Define Custom Structure";
    BOOST_TEST(dialog->windowTitle() == expected_title);
    dialog->addSMILES("CC");
    dialog->accept();
    BOOST_TEST(model.getDrawTool() == DrawTool::CUSTOM_MONOMER);
    const auto monomer = model.getValue(ModelKey::CUSTOM_MONOMER)
                             .value<std::pair<QString, ChainType>>();
    BOOST_TEST(monomer.first == "CC");
    BOOST_CHECK(monomer.second == chain_type);
    BOOST_TEST(button->isChecked());
}

/**
 * Only hovered or checked custom monomer buttons request a native menu arrow.
 */
BOOST_AUTO_TEST_CASE(custom_monomer_menu_indicator)
{
    class RecordingStyle : public QProxyStyle
    {
      public:
        mutable bool m_painted = false;
        mutable bool m_has_menu = false;

        void drawComplexControl(ComplexControl control,
                                const QStyleOptionComplex* option,
                                QPainter* painter,
                                const QWidget* widget) const override
        {
            if (control == CC_ToolButton) {
                const auto* button_option =
                    qstyleoption_cast<const QStyleOptionToolButton*>(option);
                m_painted = true;
                m_has_menu =
                    button_option->features & QStyleOptionToolButton::HasMenu;
            }
            QProxyStyle::drawComplexControl(control, option, painter, widget);
        }
    } style;
    CustomMonomerButton button;
    button.setStyle(&style);
    button.resize(100, 32);
    for (const bool enabled : {false, true}) {
        for (const bool checked : {false, true}) {
            for (const bool hovered : {false, true}) {
                button.setEnabled(enabled);
                button.setChecked(checked);
                button.setAttribute(Qt::WA_UnderMouse, hovered);
                style.m_painted = false;
                QPixmap image(button.size());
                button.render(&image);
                BOOST_TEST(style.m_painted);
                BOOST_TEST(style.m_has_menu ==
                           (checked || (enabled && hovered)));
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
    const auto core_chem_count =
        db.getMonomersByPolymerType(rdkit_extensions::ChainType::CHEM).size();
    BOOST_TEST(aa_button->isHidden() == (core_aa_count + core_chem_count == 0));
    BOOST_TEST(na_button->isHidden() == (core_chem_count == 0));

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
        BOOST_REQUIRE_EQUAL(packets.size(), (amino_acid ? core_aa_count : 0) +
                                                core_chem_count + 2);
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

        QEnterEvent enter(button->rect().center(), button->rect().center(),
                          button->mapToGlobal(button->rect().center()));
        QApplication::sendEvent(button, &enter);
        QTest::mousePress(button, Qt::LeftButton);
        BOOST_TEST(
            QTest::qWaitFor([popup]() { return popup->isVisible(); }, 2000));
        QTest::mouseRelease(button, Qt::LeftButton);
        popup->close();
        BOOST_TEST(!button->underMouse());
        BOOST_TEST(!button->isDown());
        BOOST_TEST(!button->isChecked());
        BOOST_CHECK(model.getValue(key) == previous_tool);

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
        // the page's standard tool still names the last standard button used
        if (amino_acid) {
            BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
        } else {
            BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::A);
        }
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
        const auto monomer = model.getMonomerDbMonomer();
        BOOST_TEST(monomer.symbol == symbol);
        BOOST_CHECK(monomer.chain_type ==
                    (amino_acid ? rdkit_extensions::ChainType::PEPTIDE
                                : rdkit_extensions::ChainType::RNA));
    }

    // Each page stores its own database monomer. Switching pages restores each
    // page's remembered monomer without opening its popup.
    for (const bool amino_acid : {true, false, true, false}) {
        widget
            .findChild<QAbstractButton*>(amino_acid ? "amino_monomer_btn"
                                                    : "nucleic_monomer_btn")
            ->click();
        auto* button = amino_acid ? aa_button : na_button;
        BOOST_TEST(button->isChecked());
        BOOST_TEST(!button->getPopupWidget()->isVisible());
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
        BOOST_TEST(model.getMonomerDbMonomer().symbol ==
                   (amino_acid ? "testX2" : "testN2"));
    }
    BOOST_TEST(model.getMonomerDbMonomer(MonomerToolType::AMINO_ACID).symbol ==
               "testX2");
    BOOST_TEST(
        model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID).symbol ==
        "testN2");
    BOOST_TEST(model.getValueString(ModelKey::AMINO_ACID_SYMBOL) == "A");
    BOOST_CHECK(model.getValue(ModelKey::AMINO_ACID_SYMBOL).metaType() ==
                QMetaType::fromType<QString>());
    BOOST_TEST(model.getValue(ModelKey::NUCLEIC_ACID_SYMBOL)
                   .value<NucleicAcidMutation>()
                   .symbol == "A");

    QPointer<QWidget> old_aa_popup = aa_button->getPopupWidget();
    QPointer<QWidget> old_na_popup = na_button->getPopupWidget();
    db.loadMonomersFromJson("[]");
    BOOST_TEST(old_aa_popup.isNull());
    BOOST_TEST(old_na_popup.isNull());
    BOOST_TEST(aa_button->isHidden() == (core_aa_count + core_chem_count == 0));
    BOOST_TEST(na_button->isHidden() == (core_chem_count == 0));
    BOOST_TEST(aa_button->getEnumItem() == -1);
    BOOST_TEST(na_button->getEnumItem() == -1);
    BOOST_TEST(aa_button->text() == "Unclassified");
    BOOST_TEST(na_button->text() == "Unclassified");
    // The removed monomers were cleared from both pages
    BOOST_TEST(
        model.getMonomerDbMonomer(MonomerToolType::AMINO_ACID).symbol.empty());
    BOOST_TEST(model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID)
                   .symbol.empty());
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
}

/**
 * CHEM entries are selectable from both popups and retain their identity.
 */
BOOST_AUTO_TEST_CASE(unclassified_chem_monomers)
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
    auto scene = TestScene::getScene();
    auto* model = scene->m_sketcher_model;
    model->setValue(ModelKey::TOOL_SET, ToolSet::MONOMERIC);
    MonomerToolWidget widget;
    widget.setModel(model);
    widget.show();
    const auto core_chem_count =
        db.getMonomersByPolymerType(rdkit_extensions::ChainType::CHEM).size();
    const auto core_aa_count =
        db.getMonomersByNaturalAnalog(rdkit_extensions::ChainType::PEPTIDE)["X"]
            .size();

    // Include an empty natural analog, a self analog, and a symbol matching
    // the RNA unknown slot. None may be excluded from the CHEM list.
    const auto result = db.loadMonomersFromJson(R"([
        {"symbol":"testChem","polymer_type":"CHEM","natural_analog":"",
         "smiles":"CC","name":"Chemical monomer","monomer_type":"undefined","author":"test"},
        {"symbol":"testSelf","polymer_type":"CHEM","natural_analog":"testSelf",
         "smiles":"CN","name":"Self analog","monomer_type":"undefined","author":"test"},
        {"symbol":"N","polymer_type":"CHEM","natural_analog":"N",
         "smiles":"CCC","name":"Chemical N","monomer_type":"undefined","author":"test"},
        {"symbol":"testChem","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CCN","name":"Peptide with shared symbol","monomer_type":"backbone","author":"test"},
        {"symbol":"testChem","polymer_type":"RNA","natural_analog":"N",
         "smiles":"CCO","name":"Base with shared symbol","monomer_type":"branch","author":"test"}
    ])");
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), 5);

    for (const bool amino_acid : {true, false}) {
        widget
            .findChild<QAbstractButton*>(amino_acid ? "amino_monomer_btn"
                                                    : "nucleic_monomer_btn")
            ->click();
        auto* button = widget.findChild<ModularToolButton*>(
            amino_acid ? "aa_unclassified_btn" : "na_unclassified_btn");
        BOOST_REQUIRE(button != nullptr);
        BOOST_TEST(!button->isHidden());
        auto* popup = dynamic_cast<ModularPopup*>(button->getPopupWidget());
        BOOST_REQUIRE(popup != nullptr);
        BOOST_TEST(popup->getButtonPackets().size() ==
                   core_chem_count + (amino_acid ? core_aa_count : 0) + 4);
        button->click();
        auto buttons = popup->findChildren<QToolButton*>();
        auto find_entry = [&](const QString& name) {
            for (auto* entry : buttons) {
                if (entry->toolTip() == name) {
                    return entry;
                }
            }
            return static_cast<QToolButton*>(nullptr);
        };
        auto* chem_button = find_entry("Chemical monomer");
        BOOST_REQUIRE(chem_button != nullptr);
        BOOST_REQUIRE(find_entry("Self analog") != nullptr);
        BOOST_REQUIRE(find_entry("Chemical N") != nullptr);
        chem_button->click();
        const int chem_id = button->getEnumItem();
        const auto selection = model->getMonomerDbMonomer();
        BOOST_TEST(selection.symbol == "testChem");
        BOOST_CHECK(selection.chain_type == rdkit_extensions::ChainType::CHEM);
        BOOST_CHECK(model->getDrawTool() == DrawTool::MONOMER_DB_MONOMER);

        // Switching away and back reuses the CHEM entry rather than the
        // biological entry with the same symbol.
        widget.findChild<QAbstractButton*>(amino_acid ? "ala_btn" : "na_a_btn")
            ->click();
        button->click();
        BOOST_TEST(button->getEnumItem() == chem_id);

        QGraphicsSceneMouseEvent press(QEvent::GraphicsSceneMousePress);
        QGraphicsSceneMouseEvent release(QEvent::GraphicsSceneMouseRelease);
        for (auto* event : {&press, &release}) {
            event->setScenePos(QPointF(amino_acid ? 0 : 200, 0));
            event->setButton(Qt::LeftButton);
            event->setButtons(Qt::LeftButton);
        }
        const auto old_count = scene->m_mol_model->getMol()->getNumAtoms();
        scene->mousePressEvent(&press);
        scene->mouseReleaseEvent(&release);
        const auto* mol = scene->m_mol_model->getMol();
        BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), old_count + 1);
        BOOST_CHECK(get_monomer_type(mol->getAtomWithIdx(old_count)) ==
                    MonomerType::CHEM);
        BOOST_TEST(get_monomer_res_name(mol->getAtomWithIdx(old_count)) ==
                   "testChem");
    }

    auto* aa_button =
        widget.findChild<ModularToolButton*>("aa_unclassified_btn");
    auto* na_button =
        widget.findChild<ModularToolButton*>("na_unclassified_btn");
    QPointer<QWidget> old_aa_popup = aa_button->getPopupWidget();
    QPointer<QWidget> old_na_popup = na_button->getPopupWidget();
    db.loadMonomersFromJson("[]");
    BOOST_TEST(old_aa_popup.isNull());
    BOOST_TEST(old_na_popup.isNull());
    BOOST_TEST(aa_button->isHidden() == (core_aa_count + core_chem_count == 0));
    BOOST_TEST(na_button->isHidden() == (core_chem_count == 0));
}

/**
 * Load a custom monomer database with a PEPTIDE, an RNA and a CHEM monomer
 */
static void load_mixed_db_monomers()
{
    auto result =
        rdkit_extensions::MonomerDatabase::instance().loadMonomersFromJson(R"([
        {"symbol":"pepX","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CC","name":"Peptide","monomer_type":"backbone",
         "author":"test"},
        {"symbol":"rnaN","polymer_type":"RNA","natural_analog":"N",
         "smiles":"CCC","name":"Base","monomer_type":"branch","author":"test"},
        {"symbol":"chem","polymer_type":"CHEM","natural_analog":"",
         "smiles":"CN","name":"Chem","monomer_type":"undefined",
         "author":"test"}
    ])");
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), 3);
}

/**
 * Confirm that the widget shows the given page, with its page button checked
 */
static void check_visible_page(MonomerToolWidget& widget, bool amino_acid)
{
    auto* stack = widget.findChild<QStackedWidget*>("amino_or_nucleic_stack");
    BOOST_REQUIRE(stack != nullptr);
    BOOST_TEST(stack->currentWidget()->objectName().toStdString() ==
               (amino_acid ? "amino_page" : "nucleic_page"));
    BOOST_TEST(
        widget.findChild<QAbstractButton*>("amino_monomer_btn")->isChecked() ==
        amino_acid);
    BOOST_TEST(widget.findChild<QAbstractButton*>("nucleic_monomer_btn")
                   ->isChecked() == !amino_acid);
}

/**
 * Without a canvas selection, activating a database monomer switches to the
 * monomer's page (CHEM monomers stay on the current page), activates the
 * monomer, and remembers it only on that page
 */
BOOST_AUTO_TEST_CASE(activate_monomer_db_monomer)
{
    using rdkit_extensions::MonomerID;
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    load_mixed_db_monomers();
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    const MonomerID pep{"pepX", ChainType::PEPTIDE};
    const MonomerID rna{"rnaN", ChainType::RNA};
    const MonomerID chem{"chem", ChainType::CHEM};
    const MonomerID none{"", ChainType::PEPTIDE};
    auto* aa_button =
        widget.findChild<ModularToolButton*>("aa_unclassified_btn");
    auto* na_button =
        widget.findChild<ModularToolButton*>("na_unclassified_btn");

    struct Expected {
        MonomerID monomer;
        bool amino_acid_page;
        MonomerID aa_slot;
        MonomerID na_slot;
    };
    for (const auto& expected : {
             // an RNA monomer switches from the default AA page to the NA page
             Expected{rna, false, none, rna},
             // CHEM stays on the current page
             Expected{chem, false, none, chem},
             // a peptide switches to the AA page and leaves the NA page's
             // last-used database monomer alone
             Expected{pep, true, pep, chem},
             Expected{chem, true, chem, chem},
         }) {
        BOOST_TEST_CONTEXT(expected.monomer.symbol)
        {
            widget.activateMonomerDbMonomer(expected.monomer);
            const auto page = expected.amino_acid_page
                                  ? MonomerToolType::AMINO_ACID
                                  : MonomerToolType::NUCLEIC_ACID;
            BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
            BOOST_CHECK(model.getMonomerToolType() == page);
            BOOST_CHECK(model.getMonomerDbMonomer() == expected.monomer);
            BOOST_CHECK(model.getMonomerDbMonomer(
                            MonomerToolType::AMINO_ACID) == expected.aa_slot);
            BOOST_CHECK(model.getMonomerDbMonomer(
                            MonomerToolType::NUCLEIC_ACID) == expected.na_slot);
            BOOST_TEST((model.getDrawToolToRestore(page) ==
                        DrawTool::MONOMER_DB_MONOMER));
            check_visible_page(widget, expected.amino_acid_page);
            BOOST_TEST((expected.amino_acid_page ? aa_button : na_button)
                           ->isChecked());
        }
    }
    // activating database monomers leaves the pages' standard tools alone
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::RNA_NUCLEOTIDE);

    // an empty monomer is ignored
    QSignalSpy pinged(&model, &SketcherModel::valuePinged);
    widget.activateMonomerDbMonomer(none);
    BOOST_TEST(pinged.empty());
    BOOST_CHECK(model.getMonomerDbMonomer() == chem);
}

/**
 * Each page restores its own last monomeric tool: its database monomer, or its
 * standard monomer. A page never restores a database tool without a database
 * monomer.
 */
BOOST_AUTO_TEST_CASE(monomer_db_page_switching)
{
    using rdkit_extensions::MonomerID;
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    load_mixed_db_monomers();
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    const MonomerID pep{"pepX", ChainType::PEPTIDE};
    const MonomerID rna{"rnaN", ChainType::RNA};
    auto* aa_page_button =
        widget.findChild<QAbstractButton*>("amino_monomer_btn");
    auto* na_page_button =
        widget.findChild<QAbstractButton*>("nucleic_monomer_btn");
    widget.activateMonomerDbMonomer(rna);
    widget.activateMonomerDbMonomer(pep);

    for (const bool amino_acid : {false, true, false}) {
        (amino_acid ? aa_page_button : na_page_button)->click();
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
        BOOST_CHECK(model.getMonomerDbMonomer() == (amino_acid ? pep : rna));
        check_visible_page(widget, amino_acid);
        BOOST_TEST(widget
                       .findChild<QAbstractButton*>(amino_acid
                                                        ? "aa_unclassified_btn"
                                                        : "na_unclassified_btn")
                       ->isChecked());
    }

    // a standard tool replaces the NA page's DB tool as its last monomeric
    // tool, but doesn't clear its last-used DB monomer
    auto* na_a_button = widget.findChild<QAbstractButton*>("na_a_btn");
    na_a_button->click();
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::A);
    aa_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
    BOOST_CHECK(model.getMonomerDbMonomer() == pep);
    na_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::A);
    BOOST_CHECK(model.getMonomerDbMonomer() == rna);
    BOOST_TEST(na_a_button->isChecked());

    // a page whose last tool was its database monomer tool never restores it
    // without a last-used database monomer, and restores its standard tool
    // instead
    model.setValues({{ModelKey::AMINO_ACID_DB_MONOMER,
                      QVariant::fromValue(MonomerID{"", ChainType::PEPTIDE})}});
    aa_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getMonomerToolType() == MonomerToolType::AMINO_ACID);
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
    check_visible_page(widget, true);
    BOOST_TEST(widget.findChild<QAbstractButton*>("ala_btn")->isChecked());
}

/**
 * Each page restores its own CHEM database monomer. A CHEM monomer belongs on
 * either page, so restoring it on a page switch must store it on the page being
 * switched to, not on the page being left.
 */
BOOST_AUTO_TEST_CASE(monomer_db_page_switching_chem)
{
    using rdkit_extensions::MonomerID;
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    auto result =
        rdkit_extensions::MonomerDatabase::instance().loadMonomersFromJson(R"([
        {"symbol":"CHM1","polymer_type":"CHEM","natural_analog":"",
         "smiles":"CN","name":"Chem 1","monomer_type":"undefined",
         "author":"test"},
        {"symbol":"CHM2","polymer_type":"CHEM","natural_analog":"",
         "smiles":"CO","name":"Chem 2","monomer_type":"undefined",
         "author":"test"}
    ])");
    BOOST_REQUIRE(result.second.empty());
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    const MonomerID chem1{"CHM1", ChainType::CHEM};
    const MonomerID chem2{"CHM2", ChainType::CHEM};
    auto* aa_page_button =
        widget.findChild<QAbstractButton*>("amino_monomer_btn");
    auto* na_page_button =
        widget.findChild<QAbstractButton*>("nucleic_monomer_btn");
    aa_page_button->click();
    widget.activateMonomerDbMonomer(chem2);
    na_page_button->click();
    widget.activateMonomerDbMonomer(chem1);
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::AMINO_ACID) ==
                chem2);
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID) ==
                chem1);

    for (const bool amino_acid : {true, false, true}) {
        (amino_acid ? aa_page_button : na_page_button)->click();
        BOOST_CHECK(model.getMonomerToolType() ==
                    (amino_acid ? MonomerToolType::AMINO_ACID
                                : MonomerToolType::NUCLEIC_ACID));
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
        BOOST_CHECK(model.getMonomerDbMonomer() ==
                    (amino_acid ? chem2 : chem1));
        // the other page's last-used database monomer is untouched
        BOOST_CHECK(model.getMonomerDbMonomer(
                        amino_acid ? MonomerToolType::NUCLEIC_ACID
                                   : MonomerToolType::AMINO_ACID) ==
                    (amino_acid ? chem1 : chem2));
        check_visible_page(widget, amino_acid);
    }
}

/**
 * Sketch a custom monomer through the custom monomer button's menu, as a user
 * would. The menu offers a monomer of the current page's chain type, or a CHEM
 * structure.
 */
static void sketch_custom_monomer(MonomerToolWidget& widget, bool chem,
                                  const QString& smiles)
{
    auto* button = widget.findChild<CustomMonomerButton*>("custom_monomer_btn");
    BOOST_REQUIRE(button != nullptr);
    auto* menu = button->findChild<QMenu*>();
    BOOST_REQUIRE(menu != nullptr);
    menu->actions()[chem ? 2 : 0]->trigger();
    // Earlier dialogs are deleted later, so find the one that's open
    CustomMonomerDialog* dialog = nullptr;
    for (auto* child : button->findChildren<CustomMonomerDialog*>()) {
        if (child->isVisible()) {
            dialog = child;
        }
    }
    BOOST_REQUIRE(dialog != nullptr);
    dialog->addSMILES(smiles.toStdString());
    dialog->accept();
}

/**
 * Switching to the other page and straight back restores the connection tool
 * that was active, on either page
 */
BOOST_DATA_TEST_CASE(page_round_trip_restores_connection_tool,
                     bdata::make({true, false}) *
                         bdata::make({"covalent_or_disulfide_btn",
                                      "hbond_btn"}),
                     start_on_amino_acid_page, connection_button_name)
{
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    auto* aa_page_button =
        widget.findChild<QAbstractButton*>("amino_monomer_btn");
    auto* na_page_button =
        widget.findChild<QAbstractButton*>("nucleic_monomer_btn");
    auto* start_page_button =
        start_on_amino_acid_page ? aa_page_button : na_page_button;
    auto* other_page_button =
        start_on_amino_acid_page ? na_page_button : aa_page_button;
    auto* connection_button =
        widget.findChild<QAbstractButton*>(connection_button_name);
    BOOST_REQUIRE(connection_button != nullptr);
    start_page_button->click();
    connection_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMERIC_CONNECTION);
    const auto connection_tool =
        model.getValueInt(ModelKey::MONOMERIC_CONNECTION_TOOL);

    // the other page restores its own standard tool
    other_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_TEST(!connection_button->isChecked());

    start_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMERIC_CONNECTION);
    BOOST_TEST(model.getValueInt(ModelKey::MONOMERIC_CONNECTION_TOOL) ==
               connection_tool);
    check_visible_page(widget, start_on_amino_acid_page);
    BOOST_TEST(connection_button->isChecked());
}

/**
 * Switching pages while a shared tool (select, move-rotate or erase) is active
 * keeps that tool active and doesn't disturb either page's last monomeric tool
 */
BOOST_DATA_TEST_CASE(page_round_trip_keeps_shared_tool,
                     bdata::make({DrawTool::SELECT, DrawTool::MOVE_ROTATE,
                                  DrawTool::ERASE}),
                     shared_tool)
{
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    auto* aa_page_button =
        widget.findChild<QAbstractButton*>("amino_monomer_btn");
    auto* na_page_button =
        widget.findChild<QAbstractButton*>("nucleic_monomer_btn");
    auto* hbond_button = widget.findChild<QAbstractButton*>("hbond_btn");
    auto* na_c_button = widget.findChild<QAbstractButton*>("na_c_btn");
    // the amino acid page's last tool is the H-bond tool, and the nucleic acid
    // page's is C
    model.setValue(ModelKey::DRAW_TOOL, DrawTool::MONOMER);
    hbond_button->click();
    na_page_button->click();
    na_c_button->click();
    BOOST_REQUIRE(model.getDrawTool() == DrawTool::MONOMER);

    model.setValue(ModelKey::DRAW_TOOL, shared_tool);
    for (const bool amino_acid : {true, false, true}) {
        (amino_acid ? aa_page_button : na_page_button)->click();
        BOOST_CHECK(model.getDrawTool() == shared_tool);
        BOOST_CHECK(model.getMonomerToolType() ==
                    (amino_acid ? MonomerToolType::AMINO_ACID
                                : MonomerToolType::NUCLEIC_ACID));
        check_visible_page(widget, amino_acid);
        BOOST_TEST(!hbond_button->isChecked());
        BOOST_TEST(!na_c_button->isChecked());
    }
    BOOST_CHECK(model.getDrawToolToRestore(MonomerToolType::AMINO_ACID) ==
                DrawTool::MONOMERIC_CONNECTION);
    BOOST_CHECK(model.getDrawToolToRestore(MonomerToolType::NUCLEIC_ACID) ==
                DrawTool::MONOMER);

    // a monomeric tool chosen afterwards is recorded as usual, and the other
    // page still restores its own tool
    widget.findChild<QAbstractButton*>("gly_btn")->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::GLY);
    na_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::C);
    BOOST_TEST(na_c_button->isChecked());
    aa_page_button->click();
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::GLY);
    BOOST_TEST(widget.findChild<QAbstractButton*>("gly_btn")->isChecked());
}

/**
 * Switching to the other page and straight back restores a custom monomer
 * tool, as long as the shared stored custom monomer still belongs on the page
 */
BOOST_AUTO_TEST_CASE(page_round_trip_restores_custom_monomer_tool)
{
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    auto* aa_page_button =
        widget.findChild<QAbstractButton*>("amino_monomer_btn");
    auto* na_page_button =
        widget.findChild<QAbstractButton*>("nucleic_monomer_btn");
    auto* custom_button =
        widget.findChild<QAbstractButton*>("custom_monomer_btn");
    auto* ala_button = widget.findChild<QAbstractButton*>("ala_btn");
    auto check_custom_monomer = [&](const QString& smiles,
                                    ChainType chain_type) {
        BOOST_CHECK(model.getDrawTool() == DrawTool::CUSTOM_MONOMER);
        const auto monomer = model.getCustomMonomer();
        BOOST_TEST(monomer.first == smiles);
        BOOST_CHECK(monomer.second == chain_type);
        BOOST_TEST(custom_button->isChecked());
    };

    // a custom peptide on the amino acid page
    aa_page_button->click();
    sketch_custom_monomer(widget, /*chem=*/false, "CC");
    check_custom_monomer("CC", ChainType::PEPTIDE);
    na_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_TEST(!custom_button->isChecked());
    aa_page_button->click();
    check_custom_monomer("CC", ChainType::PEPTIDE);
    check_visible_page(widget, true);

    // a custom base on the nucleic acid page replaces the stored peptide, so
    // the amino acid page falls back to its standard tool
    na_page_button->click();
    sketch_custom_monomer(widget, /*chem=*/false, "CCC");
    check_custom_monomer("CCC", ChainType::RNA);
    aa_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
    BOOST_TEST(ala_button->isChecked());
    na_page_button->click();
    check_custom_monomer("CCC", ChainType::RNA);
    check_visible_page(widget, false);

    // a CHEM custom monomer can be restored on either page
    sketch_custom_monomer(widget, /*chem=*/true, "CO");
    check_custom_monomer("CO", ChainType::CHEM);
    aa_page_button->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    na_page_button->click();
    check_custom_monomer("CO", ChainType::CHEM);
    check_visible_page(widget, false);
}

/**
 * The page follows MONOMER_TOOL_TYPE under every draw tool
 */
BOOST_DATA_TEST_CASE(monomer_page_follows_model,
                     bdata::make({DrawTool::SELECT, DrawTool::MOVE_ROTATE,
                                  DrawTool::MONOMERIC_CONNECTION,
                                  DrawTool::MONOMER}),
                     draw_tool)
{
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    model.setValue(ModelKey::DRAW_TOOL, draw_tool);
    for (const bool amino_acid : {false, true}) {
        model.setValue(ModelKey::MONOMER_TOOL_TYPE,
                       amino_acid ? MonomerToolType::AMINO_ACID
                                  : MonomerToolType::NUCLEIC_ACID);
        check_visible_page(widget, amino_acid);
    }
}

/**
 * Removing an inactive page's database monomer means that returning to the
 * page restores its standard fallback tool
 */
BOOST_AUTO_TEST_CASE(monomer_db_repair_inactive_page_fallback)
{
    using rdkit_extensions::MonomerID;
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    struct ResetDatabase {
        ~ResetDatabase()
        {
            rdkit_extensions::MonomerDatabase::instance()
                .resetMonomerDefinitions();
        }
    } reset_database;
    load_mixed_db_monomers();
    SketcherModel model;
    MonomerToolWidget widget;
    widget.setModel(&model);
    widget.show();
    widget.activateMonomerDbMonomer({"rnaN", ChainType::RNA});
    widget.activateMonomerDbMonomer({"ac", ChainType::PEPTIDE});
    auto* na_button =
        widget.findChild<ModularToolButton*>("na_unclassified_btn");
    BOOST_TEST(na_button->getEnumItem() != -1);
    db.resetMonomerDefinitions();
    // the AA page's core monomer is still available
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
    BOOST_CHECK(model.getMonomerDbMonomer() ==
                (MonomerID{"ac", ChainType::PEPTIDE}));
    BOOST_TEST(na_button->getEnumItem() == -1);
    widget.findChild<QAbstractButton*>("nucleic_monomer_btn")->click();
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::RNA_NUCLEOTIDE);
    BOOST_TEST(widget.findChild<QAbstractButton*>("na_rna_btn")->isChecked());
}

} // namespace sketcher
} // namespace schrodinger
