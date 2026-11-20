#define BOOST_TEST_MODULE test_monomer_database_popup
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QGraphicsSceneMouseEvent>
#include <QHeaderView>
#include <QHoverEvent>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QPainter>
#include <QProxyStyle>
#include <QSignalSpy>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>
#include <QToolButton>
#include <QUndoStack>

#include "../monomer_database_test_helpers.h"
#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/dialog/custom_title_bar.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.h"
#include "schrodinger/sketcher/widget/monomer_database_filter_popup.h"
#include "schrodinger/sketcher/widget/monomer_database_last_used_button.h"
#include "schrodinger/sketcher/widget/monomer_database_popup.h"
#include "schrodinger/sketcher/widget/monomer_database_table_view.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"
#include "schrodinger/sketcher/widget/multi_select_combo_box.h"
#include "schrodinger/sketcher/widget/sketcher_top_bar.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

std::ostream& operator<<(std::ostream& os, const QPoint& point)
{
    return os << "QPoint(" << point.x() << ", " << point.y() << ")";
}

std::ostream& operator<<(std::ostream& os, const QSize& size)
{
    return os << "QSize(" << size.width() << ", " << size.height() << ")";
}

namespace schrodinger
{
namespace sketcher
{

using rdkit_extensions::ChainType;
using rdkit_extensions::MonomerDatabase;
using rdkit_extensions::MonomerID;

/**
 * Custom monomers of every polymer type, including a custom override of the
 * core "am" monomer
 */
static const char* const TEST_MONOMERS_JSON = R"([
    {"symbol":"pepA","polymer_type":"PEPTIDE","natural_analog":"A",
     "smiles":"CC","name":"Peptide alanine analog","monomer_type":"backbone",
     "author":"test"},
    {"symbol":"pepX","polymer_type":"PEPTIDE","natural_analog":"X",
     "smiles":"CN","name":"Peptide without analog","monomer_type":"backbone",
     "author":"test"},
    {"symbol":"Zeta","polymer_type":"PEPTIDE","natural_analog":"G",
     "smiles":"CO","name":"Hydroxy zeta","monomer_type":"backbone",
     "author":"test"},
    {"symbol":"am","polymer_type":"PEPTIDE","natural_analog":"X",
     "smiles":"CCC","name":"Custom amide","monomer_type":"terminal",
     "author":"test"},
    {"symbol":"rnaA","polymer_type":"RNA","natural_analog":"A",
     "smiles":"CCN","name":"Base adenine analog","monomer_type":"branch",
     "author":"test"},
    {"symbol":"chem1","polymer_type":"CHEM","natural_analog":"",
     "smiles":"CCO","name":"Hydroxy structure","monomer_type":"undefined",
     "author":"test"}
])";
static constexpr std::size_t NUM_TEST_MONOMERS = 6;

/**
 * A monomer tool widget with its own model. The monomer database is reset when
 * the fixture is destroyed.
 */
struct PopupFixture {
    PopupFixture()
    {
        widget.setModel(&model);
        widget.show();
    }

    MonomerDatabasePopup* openPopup()
    {
        return open_monomer_db_popup(widget);
    }

    MonomerDatabaseLastUsedButton* lastUsedButton()
    {
        auto* button = widget.findChild<MonomerDatabaseLastUsedButton*>(
            "monomer_db_last_used_btn");
        BOOST_REQUIRE(button != nullptr);
        return button;
    }

    void clickPageButton(bool amino_acid)
    {
        widget
            .findChild<QAbstractButton*>(amino_acid ? "amino_monomer_btn"
                                                    : "nucleic_monomer_btn")
            ->click();
    }

    ResetMonomerDatabase reset_database;
    SketcherModel model;
    MonomerToolWidget widget;
};

template <typename T>
static T* find_popup_child(MonomerDatabasePopup& popup, const char* name)
{
    auto* child = popup.findChild<T*>(name);
    BOOST_REQUIRE_MESSAGE(child != nullptr, name);
    return child;
}

static QLineEdit* get_search_box(MonomerDatabasePopup& popup)
{
    return find_popup_child<QLineEdit>(popup, "monomer_db_search_le");
}

static QString get_match_count_text(MonomerDatabasePopup& popup)
{
    return find_popup_child<QLabel>(popup, "monomer_db_match_count_lbl")
        ->text();
}

static QToolButton* get_load_link(MonomerDatabasePopup& popup)
{
    return find_popup_child<QToolButton>(popup, "monomer_db_load_link");
}

/**
 * @return The table's selected rows
 */
static std::vector<int> get_selected_rows(QTableView& table)
{
    std::vector<int> rows;
    for (const auto& index : table.selectionModel()->selectedRows()) {
        rows.push_back(index.row());
    }
    return rows;
}

/**
 * @return The symbol of the table's current row, or "" if there isn't one
 */
static std::string get_current_symbol(QTableView& table)
{
    auto index = table.currentIndex();
    if (!index.isValid()) {
        return "";
    }
    return table.model()->index(index.row(), 0).data().toString().toStdString();
}

/**
 * Count the monomers that the popup chooses
 */
struct ChosenMonomers {
    ChosenMonomers(MonomerDatabasePopup& popup)
    {
        QObject::connect(
            &popup, &MonomerDatabasePopup::monomerChosen,
            [this](const MonomerID& monomer) { monomers.push_back(monomer); });
    }
    std::vector<MonomerID> monomers;
};

/**
 * With no custom database, the table shows the core monomers with no natural
 * analog and the link offers to load a database. Once custom monomers are
 * loaded, it offers to update the database instead.
 */
BOOST_AUTO_TEST_CASE(load_link_text)
{
    PopupFixture fixture;
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"ac", "am"}),
               boost::test_tools::per_element());
    BOOST_TEST(get_match_count_text(*popup) == "2 matches");
    auto* load_link = get_load_link(*popup);
    BOOST_TEST(load_link->isVisible());
    BOOST_TEST(load_link->text() == "Load Database...");

    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    BOOST_TEST(load_link->text() == "Update Database...");
    MonomerDatabase::instance().resetMonomerDefinitions();
    BOOST_TEST(load_link->text() == "Load Database...");
    // an empty custom database has no custom monomers
    load_custom_monomers("[]", 0);
    BOOST_TEST(load_link->text() == "Load Database...");
}

/**
 * Hiding the load link applies whether or not the popup exists yet, and
 * SketcherWidget::setLoadMonomerDatabaseVisible() hides both the link and the
 * import menu's action
 */
BOOST_AUTO_TEST_CASE(load_link_visibility)
{
    {
        PopupFixture fixture;
        fixture.widget.setLoadMonomerDatabaseVisible(false);
        auto* popup = fixture.openPopup();
        BOOST_TEST(!get_load_link(*popup)->isVisible());
        fixture.widget.setLoadMonomerDatabaseVisible(true);
        BOOST_TEST(get_load_link(*popup)->isVisible());
        fixture.widget.setLoadMonomerDatabaseVisible(false);
        BOOST_TEST(!get_load_link(*popup)->isVisible());
    }

    TestSketcherWidget sk;
    sk.setInterfaceType(InterfaceType::ATOMISTIC_OR_MONOMERIC);
    sk.m_sketcher_model->setValue(ModelKey::TOOL_SET, ToolSet::MONOMERIC);
    sk.setLoadMonomerDatabaseVisible(false);
    auto* top_bar = sk.findChild<SketcherTopBar*>();
    BOOST_REQUIRE(top_bar != nullptr);
    auto* import_menu = top_bar->findChild<ImportMenu*>();
    BOOST_REQUIRE(import_menu != nullptr);
    auto* action = import_menu->m_load_monomer_database_act;
    BOOST_TEST(!action->isVisible());
    auto* popup = open_monomer_db_popup(sk);
    BOOST_TEST(!get_load_link(*popup)->isVisible());
    popup->hide();
    sk.setLoadMonomerDatabaseVisible(true);
    BOOST_TEST(action->isVisible());
    popup = open_monomer_db_popup(sk);
    BOOST_TEST(get_load_link(*popup)->isVisible());
    popup->hide();
}

/**
 * @return The visible file chooser, which must exist
 */
static QFileDialog* get_file_chooser(QWidget& parent)
{
    QFileDialog* chooser = nullptr;
    BOOST_TEST(QTest::qWaitFor(
        [&]() {
            for (auto* dialog : parent.findChildren<QFileDialog*>()) {
                if (dialog->isVisible()) {
                    chooser = dialog;
                }
            }
            return chooser != nullptr;
        },
        2000));
    BOOST_REQUIRE(chooser != nullptr);
    return chooser;
}

/**
 * Close the visible result dialog, which must have the given title
 */
static void close_result_dialog(QWidget& parent, const QString& title)
{
    MessageBoxDialog* result = nullptr;
    for (auto* dialog : parent.findChildren<MessageBoxDialog*>()) {
        if (dialog->isVisible()) {
            result = dialog;
        }
    }
    BOOST_REQUIRE(result != nullptr);
    BOOST_TEST(result->windowTitle() == title);
    result->accept();
    QApplication::processEvents();
}

/**
 * The load link closes the popup before the file chooser opens, and the popup
 * stays closed after success, cancellation, a partial import or a failure.
 * Reopening it on the same page keeps the search text and sort and shows the
 * new rows, while reopening it on the other page resets the search text.
 */
BOOST_AUTO_TEST_CASE(load_link_lifecycle)
{
    PopupFixture fixture;
    QTemporaryDir temp_dir;
    BOOST_REQUIRE(temp_dir.isValid());
    auto write_file = [&temp_dir](const QString& name, const QByteArray& text) {
        auto path = temp_dir.filePath(name);
        QFile file(path);
        BOOST_REQUIRE(file.open(QIODevice::WriteOnly));
        file.write(text);
        return path;
    };
    const auto valid_path =
        write_file("valid.json", QByteArray(TEST_MONOMERS_JSON));
    const auto partial_path = write_file("partial.json", R"([
        {"symbol":"pepB","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CC","name":"Peptide B","monomer_type":"backbone",
         "author":"test"},
        {"symbol":"broken","polymer_type":"PEPTIDE"}
    ])");
    const auto invalid_path = write_file("invalid.json", "not JSON");

    auto* popup = fixture.openPopup();
    auto* search_box = get_search_box(*popup);
    QTest::keyClicks(search_box, "a");
    auto* table = get_monomer_db_table(*popup);
    table->sortByColumn(static_cast<int>(MonomerDbColumn::NAME),
                        Qt::DescendingOrder);
    const MonomerDbSort expected_sort{MonomerDbColumn::NAME, false};
    BOOST_REQUIRE(popup->getSort() == expected_sort);

    struct Outcome {
        QString path;         // empty to cancel
        QString result_title; // empty if no result dialog is expected
    };
    for (const auto& outcome :
         {Outcome{"", ""}, Outcome{invalid_path, "Monomer Database Error"},
          Outcome{partial_path, "Monomer Database Error"},
          Outcome{valid_path, "Monomer Database"}}) {
        BOOST_TEST_CONTEXT("path=" << outcome.path)
        {
            if (!popup->isVisible()) {
                popup = fixture.openPopup();
            }
            get_load_link(*popup)->click();
            BOOST_TEST(!popup->isVisible());
            auto* chooser = get_file_chooser(fixture.widget);
            BOOST_TEST(!popup->isVisible());
            if (outcome.path.isEmpty()) {
                chooser->reject();
            } else {
                // The chooser reports the chosen file with this signal
                emit chooser->fileSelected(outcome.path);
                chooser->close();
            }
            QApplication::processEvents();
            BOOST_TEST(!popup->isVisible());
            if (!outcome.result_title.isEmpty()) {
                close_result_dialog(fixture.widget, outcome.result_title);
            }
            BOOST_TEST(!popup->isVisible());
        }
    }

    // reopening on the same page keeps the search text and sort, and shows the
    // monomers from the new database
    popup = fixture.openPopup();
    BOOST_TEST(search_box->text() == "a");
    BOOST_CHECK(popup->getSort() == expected_sort);
    BOOST_TEST(table->horizontalHeader()->sortIndicatorSection() ==
               static_cast<int>(MonomerDbColumn::NAME));
    BOOST_CHECK(table->horizontalHeader()->sortIndicatorOrder() ==
                Qt::DescendingOrder);
    BOOST_TEST(find_monomer_db_row(*table, "pepA").has_value());
    BOOST_TEST(!find_monomer_db_row(*table, "pepB").has_value());
    popup->hide();

    // changing the page while the popup is closed resets the search text
    fixture.clickPageButton(false);
    popup = fixture.openPopup();
    BOOST_TEST(search_box->text().isEmpty());
    BOOST_CHECK(popup->getSort() == expected_sort);
    popup->hide();
}

/**
 * The table holds the custom monomers plus the core monomers with no natural
 * analog, and a custom monomer replaces a core monomer with the same symbol.
 * Each page's default filter shows its own polymer type and structures.
 */
BOOST_AUTO_TEST_CASE(table_contents)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"ac", "am", "chem1", "pepA",
                                             "pepX", "Zeta"}),
               boost::test_tools::per_element());
    BOOST_TEST(!table->isColumnHidden(
        static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));
    auto am_row = find_monomer_db_row(*table, "am");
    BOOST_REQUIRE(am_row.has_value());
    BOOST_TEST(table->model()->index(*am_row, 1).data().toString() ==
               "Custom amide");
    BOOST_TEST(get_match_count_text(*popup) == "6 matches");
    popup->hide();

    fixture.clickPageButton(false);
    popup = fixture.openPopup();
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"chem1", "rnaA"}),
               boost::test_tools::per_element());
    BOOST_CHECK(
        popup->getFilter().getPolymerTypes() ==
        (MonomerDbPolymerType::NUCLEIC_ACID | MonomerDbPolymerType::CHEM));
    BOOST_TEST(!table->isColumnHidden(
        static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));
    popup->hide();
}

/**
 * The popup's title bar shows its title and a close button, and the popup's
 * outline matches the title bar's color
 */
BOOST_AUTO_TEST_CASE(title_bar)
{
    PopupFixture fixture;
    auto* popup = fixture.openPopup();
    auto* title_bar =
        find_popup_child<CustomTitleBar>(*popup, "monomer_db_title_bar");
    BOOST_TEST(title_bar->m_title_label->text() == "Monomer Database");
    auto* close_button =
        find_popup_child<QToolButton>(*popup, "monomer_db_close_btn");
    BOOST_TEST(title_bar->getCloseButton() == close_button);

    const auto image = popup->grab().toImage();
    for (const auto& point : {QPoint(0, image.height() / 2),
                              QPoint(image.width() - 1, image.height() / 2),
                              QPoint(image.width() / 2, image.height() - 1)}) {
        BOOST_TEST(image.pixelColor(point).name() == POPUP_TITLE_BAR_COLOR);
    }
    BOOST_TEST(image.pixelColor(title_bar->geometry().center()).name() ==
               POPUP_TITLE_BAR_COLOR);

    close_button->click();
    BOOST_TEST(!popup->isVisible());
}

/**
 * The search text and filter are kept when the popup reopens on the page that
 * it closed on, and reset when it opens on a different page. The page is
 * recorded when the popup closes, however it closes.
 */
BOOST_AUTO_TEST_CASE(reset_rule)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* search_box = get_search_box(*popup);
    QTest::keyClicks(search_box, "pep");
    BOOST_TEST(get_match_count_text(*popup) == "2 matches");

    // each way of closing the popup records the page
    auto close_with_x = [&]() {
        find_popup_child<QToolButton>(*popup, "monomer_db_close_btn")->click();
    };
    auto close_with_link = [&]() {
        find_popup_child<QToolButton>(*popup, "monomer_db_close_link")->click();
    };
    auto close_with_escape = [&]() {
        QTest::keyClick(search_box, Qt::Key_Escape);
    };
    auto close_with_db_button = [&]() {
        fixture.widget.findChild<QToolButton*>("monomer_db_btn")->click();
    };
    for (const auto& close : std::vector<std::function<void()>>{
             close_with_x, close_with_link, close_with_escape,
             close_with_db_button}) {
        close();
        BOOST_TEST(!popup->isVisible());
        popup = fixture.openPopup();
        BOOST_TEST(search_box->text() == "pep");
    }

    // a page change while the popup is open doesn't reset anything, since the
    // page is recorded when the popup closes
    fixture.widget.activateMonomerDbMonomer({"rnaA", ChainType::RNA});
    BOOST_TEST(popup->isVisible());
    popup->hide();
    popup = fixture.openPopup();
    BOOST_TEST(search_box->text() == "pep");
    BOOST_CHECK(popup->getFilter().getPolymerTypes() ==
                (MonomerDbPolymerType::PEPTIDE | MonomerDbPolymerType::CHEM));
    popup->hide();

    // a page change while the popup is closed resets the search text and the
    // filter
    fixture.clickPageButton(true);
    popup = fixture.openPopup();
    BOOST_TEST(search_box->text().isEmpty());
    BOOST_CHECK(popup->getFilter().getPolymerTypes() ==
                (MonomerDbPolymerType::PEPTIDE | MonomerDbPolymerType::CHEM));
    BOOST_TEST(get_match_count_text(*popup) == "6 matches");
    popup->hide();
}

/**
 * @return Whether any pixel of the image is close to the given color
 */
static bool contains_color(const QImage& image, const QColor& color)
{
    for (int y = 0; y < image.height(); ++y) {
        for (int x = 0; x < image.width(); ++x) {
            const QColor pixel = image.pixelColor(x, y);
            if (std::abs(pixel.red() - color.red()) < 40 &&
                std::abs(pixel.green() - color.green()) < 40 &&
                std::abs(pixel.blue() - color.blue()) < 40) {
                return true;
            }
        }
    }
    return false;
}

/**
 * The search narrows the rows by symbol or name, and a search with no matches
 * shows a message in the empty table
 */
BOOST_AUTO_TEST_CASE(search)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    auto* search_box = get_search_box(*popup);
    // "hydroxy" matches the names of Zeta and chem1, case-insensitively
    QTest::keyClicks(search_box, "HYDROXY");
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"chem1", "Zeta"}),
               boost::test_tools::per_element());
    BOOST_TEST(get_match_count_text(*popup) == "2 matches");
    search_box->setText("zet");
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"Zeta"}),
               boost::test_tools::per_element());
    BOOST_TEST(get_match_count_text(*popup) == "1 match");

    search_box->setText("no such monomer");
    BOOST_TEST(table->model()->rowCount() == 0);
    BOOST_TEST(get_match_count_text(*popup) == "0 matches");
    auto image = table->viewport()->grab().toImage();
    BOOST_TEST(contains_color(image, QColor(INFO_TEXT_COLOR)));
    search_box->clear();
    BOOST_TEST(table->model()->rowCount() == 6);
    BOOST_TEST(get_match_count_text(*popup) == "6 matches");
}

/**
 * The rows are first sorted by symbol, and clicking a header sorts by its
 * column. The header's indicator and the stored sort always agree.
 */
BOOST_AUTO_TEST_CASE(sorting)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    auto* header = table->horizontalHeader();
    auto check_indicator = [&](MonomerDbColumn column, bool ascending) {
        BOOST_CHECK(popup->getSort() == (MonomerDbSort{column, ascending}));
        BOOST_TEST(header->sortIndicatorSection() == static_cast<int>(column));
        BOOST_CHECK(header->sortIndicatorOrder() ==
                    (ascending ? Qt::AscendingOrder : Qt::DescendingOrder));
    };
    check_indicator(MonomerDbColumn::SYMBOL, true);
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"ac", "am", "chem1", "pepA",
                                             "pepX", "Zeta"}),
               boost::test_tools::per_element());

    auto click_header = [&](MonomerDbColumn column) {
        const int section = static_cast<int>(column);
        QTest::mouseClick(header->viewport(), Qt::LeftButton, Qt::NoModifier,
                          QPoint(header->sectionViewportPosition(section) +
                                     header->sectionSize(section) / 2,
                                 header->height() / 2));
    };
    click_header(MonomerDbColumn::NAME);
    check_indicator(MonomerDbColumn::NAME, true);
    // Custom amide, Hydroxy structure, Hydroxy zeta, Peptide alanine analog,
    // Peptide without analog, and the core ac's name
    auto symbols = get_monomer_db_symbols(*table);
    BOOST_TEST(symbols.front() != "ac");
    auto names = std::vector<std::string>();
    for (int row = 0; row < table->model()->rowCount(); ++row) {
        names.push_back(table->model()
                            ->index(row, 1)
                            .data()
                            .toString()
                            .toLower()
                            .toStdString());
    }
    BOOST_TEST(std::is_sorted(names.begin(), names.end()));
    click_header(MonomerDbColumn::NAME);
    check_indicator(MonomerDbColumn::NAME, false);
    BOOST_TEST(get_monomer_db_symbols(*table).back() == symbols.front());

    // reopening the popup keeps the sort
    popup->hide();
    popup = fixture.openPopup();
    check_indicator(MonomerDbColumn::NAME, false);
    popup->hide();
}

/**
 * Opening the popup never changes the tool. In tool mode, choosing a row with
 * the mouse or keyboard activates its monomer and keeps the popup open.
 */
BOOST_AUTO_TEST_CASE(choose_row_tool_mode)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto& model = fixture.model;
    QSignalSpy pinged(&model, &SketcherModel::valuePinged);
    const auto draw_tool = model.getDrawTool();
    auto* popup = fixture.openPopup();
    BOOST_TEST(!popup->isMutationMode());
    BOOST_TEST(pinged.empty());
    BOOST_CHECK(model.getDrawTool() == draw_tool);
    auto* table = get_monomer_db_table(*popup);
    BOOST_TEST(get_selected_rows(*table).empty());

    click_monomer_db_row(*table, "pepX");
    BOOST_TEST(popup->isVisible());
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
    BOOST_CHECK(model.getMonomerDbMonomer() ==
                (MonomerID{"pepX", ChainType::PEPTIDE}));
    auto* last_used_button = fixture.lastUsedButton();
    BOOST_TEST(last_used_button->text() == "pepX");
    BOOST_TEST(last_used_button->isChecked());
    BOOST_TEST(get_selected_rows(*table) ==
                   std::vector<int>{*find_monomer_db_row(*table, "pepX")},
               boost::test_tools::per_element());

    // keyboard navigation activates each row it moves to
    QTest::keyClick(table, Qt::Key_Down);
    BOOST_CHECK(model.getMonomerDbMonomer() ==
                (MonomerID{"Zeta", ChainType::PEPTIDE}));
    BOOST_TEST(last_used_button->text() == "Zeta");
    QTest::keyClick(table, Qt::Key_Up);
    BOOST_CHECK(model.getMonomerDbMonomer() ==
                (MonomerID{"pepX", ChainType::PEPTIDE}));
    BOOST_TEST(popup->isVisible());
    popup->hide();
}

/**
 * The table selection mirrors the model's active database monomer, without
 * choosing anything
 */
BOOST_AUTO_TEST_CASE(selection_mirroring)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto& model = fixture.model;
    auto* popup = fixture.openPopup();
    ChosenMonomers chosen(*popup);
    auto* table = get_monomer_db_table(*popup);
    auto row_of = [&](const std::string& symbol) {
        auto row = find_monomer_db_row(*table, symbol);
        BOOST_REQUIRE(row.has_value());
        return *row;
    };

    // a change made elsewhere moves the selection
    fixture.widget.activateMonomerDbMonomer({"pepA", ChainType::PEPTIDE});
    BOOST_TEST(get_selected_rows(*table) == std::vector<int>{row_of("pepA")},
               boost::test_tools::per_element());
    BOOST_TEST(get_current_symbol(*table) == "pepA");
    fixture.lastUsedButton()->click();
    BOOST_TEST(get_selected_rows(*table) == std::vector<int>{row_of("pepA")},
               boost::test_tools::per_element());

    // a non-database tool clears the selection but keeps the current index
    fixture.widget.findChild<QAbstractButton*>("ala_btn")->click();
    BOOST_TEST(get_selected_rows(*table).empty());
    BOOST_TEST(get_current_symbol(*table) == "pepA");

    // a monomer that isn't in the table selects nothing
    fixture.widget.activateMonomerDbMonomer({"A", ChainType::PEPTIDE});
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
    BOOST_TEST(get_selected_rows(*table).empty());
    // nor does one that isn't in the database at all
    model.setValues(
        {{ModelKey::AMINO_ACID_DB_MONOMER,
          QVariant::fromValue(MonomerID{"missing", ChainType::PEPTIDE})}});
    BOOST_TEST(get_selected_rows(*table).empty());

    // synchronizing never chooses a monomer
    BOOST_TEST(chosen.monomers.empty());
    popup->hide();
}

/**
 * Enter, Return and Space activate the current row in tool mode, even when a
 * standard tool replaced its monomer while the popup was open, and never close
 * the popup
 */
BOOST_AUTO_TEST_CASE(commit_keys_tool_mode)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto& model = fixture.model;
    for (const auto key : {Qt::Key_Return, Qt::Key_Enter, Qt::Key_Space}) {
        auto* popup = fixture.openPopup();
        auto* table = get_monomer_db_table(*popup);
        click_monomer_db_row(*table, "pepA");
        fixture.widget.findChild<QAbstractButton*>("ala_btn")->click();
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
        BOOST_TEST(get_selected_rows(*table).empty());
        BOOST_TEST(get_current_symbol(*table) == "pepA");
        table->setFocus();
        QTest::keyClick(table, key);
        BOOST_TEST(popup->isVisible());
        BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
        BOOST_CHECK(model.getMonomerDbMonomer() ==
                    (MonomerID{"pepA", ChainType::PEPTIDE}));
        // pressing it again on the active monomer does nothing else
        QTest::keyClick(table, key);
        BOOST_TEST(popup->isVisible());
        popup->hide();
    }
}

/**
 * The current row is highlighted when it isn't selected, and the selected row
 * uses the selection color
 */
BOOST_AUTO_TEST_CASE(current_row_highlight)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    auto row_color = [&](const std::string& symbol) {
        auto row = find_monomer_db_row(*table, symbol);
        BOOST_REQUIRE(row.has_value());
        auto image = table->viewport()->grab().toImage();
        // sample the Name column, away from the text
        auto rect = table->visualRect(table->model()->index(*row, 1));
        return image.pixelColor(rect.right() - 2, rect.center().y());
    };
    // the selected row
    click_monomer_db_row(*table, "pepA");
    BOOST_TEST(row_color("pepA").name() == CHECKED_BACKGROUND_COLOR);
    // a current row that isn't selected
    fixture.widget.findChild<QAbstractButton*>("ala_btn")->click();
    BOOST_TEST(get_current_symbol(*table) == "pepA");
    BOOST_TEST(row_color("pepA").name() == CURRENT_TABLE_ROW_COLOR);
    BOOST_TEST(row_color("pepX").name() != CURRENT_TABLE_ROW_COLOR);
    popup->hide();
}

/**
 * Reopening the popup doesn't keep the row that was current when it closed,
 * unless that row's monomer is still the active tool, which selects it
 */
BOOST_AUTO_TEST_CASE(current_row_cleared_on_reopen)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    click_monomer_db_row(*table, "pepA");
    popup->hide();

    // the database monomer is still the active tool
    popup = fixture.openPopup();
    BOOST_TEST(get_current_symbol(*table) == "pepA");
    BOOST_TEST(get_selected_rows(*table).size() == 1);
    popup->hide();

    // a standard tool replaced it
    fixture.widget.findChild<QAbstractButton*>("ala_btn")->click();
    popup = fixture.openPopup();
    BOOST_TEST(get_current_symbol(*table).empty());
    BOOST_TEST(get_selected_rows(*table).empty());
    popup->hide();
}

/**
 * Reopening the popup doesn't make a row current, even when the table had the
 * focus when the popup closed. Showing a popup gives the focus back to the
 * widget that last had it, and a table that gets the focus without a current
 * index makes its first row current.
 */
BOOST_AUTO_TEST_CASE(current_row_cleared_on_reopen_after_table_focus)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    // Closing a popup only moves the focus out of it when another window is
    // active, so a table that still has the focus wouldn't get it again
    fixture.widget.activateWindow();
    BOOST_REQUIRE(QTest::qWaitForWindowActive(&fixture.widget));
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    click_monomer_db_row(*table, "pepA");
    // A real click on a row gives the table the focus, but the synthesized one
    // above doesn't
    table->setFocus(Qt::MouseFocusReason);
    popup->hide();
    BOOST_TEST(!table->hasFocus());

    fixture.widget.findChild<QAbstractButton*>("ala_btn")->click();
    popup = fixture.openPopup();
    BOOST_TEST(get_current_symbol(*table).empty());
    BOOST_TEST(get_selected_rows(*table).empty());
    BOOST_TEST(get_search_box(*popup)->hasFocus());
    popup->hide();
}

/**
 * A style that records the state of each table cell it draws: its row
 * background and its contents
 */
class CellStateRecordingStyle : public QProxyStyle
{
  public:
    void drawPrimitive(PrimitiveElement element, const QStyleOption* option,
                       QPainter* painter, const QWidget* widget) const override
    {
        if (element == PE_PanelItemViewRow) {
            m_row_states.push_back(option->state);
        }
        QProxyStyle::drawPrimitive(element, option, painter, widget);
    }

    void drawControl(ControlElement element, const QStyleOption* option,
                     QPainter* painter, const QWidget* widget) const override
    {
        if (element == CE_ItemViewItem) {
            m_item_states.push_back(option->state);
        }
        QProxyStyle::drawControl(element, option, painter, widget);
    }

    mutable std::vector<QStyle::State> m_row_states;
    mutable std::vector<QStyle::State> m_item_states;
};

/**
 * The cell under the cursor is never drawn as hovered, so styles such as the
 * Windows style don't highlight that one cell. QTableView passes the hover
 * state to the style when it paints the cell's row background, as well as
 * through the delegate.
 */
BOOST_AUTO_TEST_CASE(no_cell_hover_highlight)
{
    PopupFixture fixture;
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    BOOST_REQUIRE(table->model()->rowCount() > 0);
    // Record what the table's own style passes on to its base style, through
    // the popup's style sheet. The table's style takes ownership of its base.
    auto* table_style = table->findChild<QProxyStyle*>();
    BOOST_REQUIRE(table_style != nullptr);
    auto* recording_style = new CellStateRecordingStyle();
    table_style->setBaseStyle(recording_style);

    // hover over a cell that isn't current or selected
    const auto cell_center =
        table->visualRect(table->model()->index(0, 1)).center();
    QHoverEvent hover(QEvent::HoverMove, cell_center,
                      table->viewport()->mapToGlobal(cell_center), cell_center);
    QApplication::sendEvent(table->viewport(), &hover);
    table->viewport()->grab();

    BOOST_TEST(!recording_style->m_row_states.empty());
    BOOST_TEST(!recording_style->m_item_states.empty());
    for (const auto& states :
         {recording_style->m_row_states, recording_style->m_item_states}) {
        for (auto state : states) {
            BOOST_TEST(!(state & QStyle::State_MouseOver));
        }
    }
    popup->hide();
}

/**
 * The current row is preserved by monomer ID across sorting, searching and
 * database refreshes while it's still in the table, and cleared once it isn't.
 * None of this chooses a monomer.
 */
BOOST_AUTO_TEST_CASE(current_index_across_resets)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    ChosenMonomers chosen(*popup);
    auto* table = get_monomer_db_table(*popup);
    auto* search_box = get_search_box(*popup);
    auto pep_x_row = *find_monomer_db_row(*table, "pepX");
    table->selectionModel()->setCurrentIndex(
        table->model()->index(pep_x_row, 0), QItemSelectionModel::NoUpdate);
    BOOST_TEST(get_current_symbol(*table) == "pepX");

    table->sortByColumn(static_cast<int>(MonomerDbColumn::SYMBOL),
                        Qt::DescendingOrder);
    BOOST_TEST(get_current_symbol(*table) == "pepX");
    BOOST_TEST(table->currentIndex().row() != pep_x_row);
    search_box->setText("pep");
    BOOST_TEST(get_current_symbol(*table) == "pepX");
    // a database refresh that keeps the monomer
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    BOOST_TEST(get_current_symbol(*table) == "pepX");

    // once the monomer is filtered out, the current index is cleared rather
    // than moved to another monomer, and it isn't restored afterwards
    search_box->setText("pepA");
    BOOST_TEST(get_current_symbol(*table) == "");
    search_box->clear();
    BOOST_TEST(get_current_symbol(*table) == "");
    BOOST_TEST(chosen.monomers.empty());
    popup->hide();
}

/**
 * A database change refreshes the table, even while the popup is hidden, and
 * removing the active database monomer blanks the last-used button
 */
BOOST_AUTO_TEST_CASE(database_reload)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    click_monomer_db_row(*table, "Zeta");
    BOOST_TEST(fixture.lastUsedButton()->text() == "Zeta");
    load_custom_monomers(R"([
        {"symbol":"newPep","polymer_type":"PEPTIDE","natural_analog":"X",
         "smiles":"CC","name":"New peptide","monomer_type":"backbone",
         "author":"test"}
    ])",
                         1);
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"ac", "am", "newPep"}),
               boost::test_tools::per_element());
    BOOST_TEST(get_selected_rows(*table).empty());
    BOOST_TEST(fixture.lastUsedButton()->text().isEmpty());
    BOOST_CHECK(fixture.model.getDrawTool() == DrawTool::MONOMER);

    popup->hide();
    MonomerDatabase::instance().resetMonomerDefinitions();
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"ac", "am"}),
               boost::test_tools::per_element());
}

/**
 * The core monomers with no natural analog can be chosen and drawn
 */
BOOST_DATA_TEST_CASE(choose_and_draw_core_monomer,
                     boost::unit_test::data::make({"ac", "am"}), symbol)
{
    ResetMonomerDatabase reset_database;
    auto scene = TestScene::getScene();
    auto* model = scene->m_sketcher_model;
    model->setValue(ModelKey::TOOL_SET, ToolSet::MONOMERIC);
    MonomerToolWidget widget;
    widget.setModel(model);
    widget.show();
    auto* popup = open_monomer_db_popup(widget);
    click_monomer_db_row(*get_monomer_db_table(*popup), symbol);
    popup->hide();
    BOOST_CHECK(model->getMonomerDbMonomer() ==
                (MonomerID{symbol, ChainType::PEPTIDE}));

    QGraphicsSceneMouseEvent press(QEvent::GraphicsSceneMousePress);
    QGraphicsSceneMouseEvent release(QEvent::GraphicsSceneMouseRelease);
    for (auto* event : {&press, &release}) {
        event->setScenePos(QPointF(0, 0));
        event->setButton(Qt::LeftButton);
        event->setButtons(Qt::LeftButton);
    }
    scene->mousePressEvent(&press);
    scene->mouseReleaseEvent(&release);
    const auto* mol = scene->m_mol_model->getMol();
    BOOST_REQUIRE_EQUAL(mol->getNumAtoms(), 1);
    BOOST_TEST(get_monomer_res_name(mol->getAtomWithIdx(0)) == symbol);
}

/**
 * A sketcher with a monomeric peptide that's selected, so that the Monomer
 * Database popup opens in mutation mode
 */
struct MutationFixture {
    MutationFixture()
    {
        load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
        sk.setInterfaceType(InterfaceType::ATOMISTIC_OR_MONOMERIC);
        sk.addFromString("PEPTIDE1{A}$$$$V2.0");
        sk.m_sketcher_model->setValue(ModelKey::TOOL_SET, ToolSet::MONOMERIC);
        sk.m_mol_model->selectAll();
        BOOST_REQUIRE(sk.m_sketcher_model->hasActiveSelection());
    }

    std::string getHelm()
    {
        return get_mol_text(sk.m_mol_model, rdkit_extensions::Format::HELM);
    }

    ResetMonomerDatabase reset_database;
    TestSketcherWidget sk;
};

/**
 * In mutation mode, keyboard navigation only moves the current index, and a
 * click, Enter, Return or Space mutates the selection once per commit. The
 * table never gains a selected row and the popup stays open.
 */
BOOST_AUTO_TEST_CASE(mutation_mode)
{
    MutationFixture fixture;
    auto& sk = fixture.sk;
    auto* sketcher_model = sk.m_sketcher_model;
    const auto draw_tool = sketcher_model->getDrawTool();
    auto* popup = open_monomer_db_popup(sk);
    BOOST_TEST(popup->isMutationMode());
    ChosenMonomers chosen(*popup);
    auto* table = get_monomer_db_table(*popup);
    auto* last_used_button = sk.findChild<MonomerDatabaseLastUsedButton*>(
        "monomer_db_last_used_btn");
    BOOST_REQUIRE(last_used_button != nullptr);
    const auto initial_undo_count = sk.m_undo_stack->count();

    // keyboard navigation doesn't mutate, add an undo entry, update the
    // last-used monomer or switch pages
    table->setFocus();
    for (int i = 0; i < 2; ++i) {
        const auto previous_symbol = get_current_symbol(*table);
        QTest::keyClick(table, Qt::Key_Down);
        BOOST_TEST(!get_current_symbol(*table).empty());
        BOOST_TEST(get_current_symbol(*table) != previous_symbol);
    }
    BOOST_TEST(fixture.getHelm() == "PEPTIDE1{A}$$$$V2.0");
    BOOST_TEST(sk.m_undo_stack->count() == initial_undo_count);
    BOOST_TEST(last_used_button->text().isEmpty());
    BOOST_CHECK(sketcher_model->getMonomerToolType() ==
                MonomerToolType::AMINO_ACID);
    BOOST_TEST(chosen.monomers.empty());
    BOOST_TEST(get_selected_rows(*table).empty());

    // a click commits once and keeps the clicked row current
    click_monomer_db_row(*table, "pepA");
    BOOST_TEST(chosen.monomers.size() == 1);
    BOOST_TEST(fixture.getHelm() == "PEPTIDE1{[pepA]}$$$$V2.0");
    BOOST_TEST(sk.m_undo_stack->count() == initial_undo_count + 1);
    BOOST_TEST(popup->isVisible());
    BOOST_TEST(get_selected_rows(*table).empty());
    BOOST_TEST(get_current_symbol(*table) == "pepA");
    BOOST_TEST(last_used_button->text() == "pepA");
    BOOST_TEST(!last_used_button->isChecked());
    BOOST_CHECK(sketcher_model->getDrawTool() == draw_tool);
    BOOST_TEST(sketcher_model->hasActiveSelection());

    // the same row can be committed again without a row change
    click_monomer_db_row(*table, "pepA");
    BOOST_TEST(chosen.monomers.size() == 2);
    BOOST_TEST(get_selected_rows(*table).empty());

    // commit keys use the current row, which is moved with the arrow keys.
    // The rows are ac, am, chem1, pepA, pepX and Zeta.
    struct Commit {
        Qt::Key move_key;
        Qt::Key commit_key;
        std::string symbol;
    };
    std::size_t expected_commits = 2;
    for (const auto& [move_key, key, symbol] :
         {Commit{Qt::Key_Down, Qt::Key_Return, "pepX"},
          Commit{Qt::Key_Down, Qt::Key_Enter, "Zeta"},
          Commit{Qt::Key_Up, Qt::Key_Space, "pepX"}}) {
        BOOST_TEST_CONTEXT("key=" << key)
        {
            table->setFocus();
            QTest::keyClick(table, move_key);
            BOOST_TEST(get_current_symbol(*table) == symbol);
            BOOST_TEST(chosen.monomers.size() == expected_commits);
            QTest::keyClick(table, key);
            ++expected_commits;
            BOOST_TEST(chosen.monomers.size() == expected_commits);
            BOOST_TEST(fixture.getHelm() ==
                       "PEPTIDE1{[" + symbol + "]}$$$$V2.0");
            BOOST_TEST(popup->isVisible());
            BOOST_TEST(get_selected_rows(*table).empty());
            BOOST_TEST(get_current_symbol(*table) == symbol);
        }
    }
    // an auto-repeated key doesn't commit again
    QKeyEvent repeat(QEvent::KeyPress, Qt::Key_Return, Qt::NoModifier, "\r",
                     /*autorep=*/true);
    QApplication::sendEvent(table, &repeat);
    BOOST_TEST(chosen.monomers.size() == expected_commits);

    // commit keys in the search box don't commit a row
    auto* search_box = get_search_box(*popup);
    for (const auto key : {Qt::Key_Return, Qt::Key_Enter, Qt::Key_Space}) {
        QTest::keyClick(search_box, key);
    }
    BOOST_TEST(chosen.monomers.size() == expected_commits);
    BOOST_TEST(popup->isVisible());
    popup->hide();
}

/**
 * The popup's mode is fixed while it's open, so a programmatic change that
 * adds or clears the canvas selection closes it. A change that keeps a
 * selection doesn't.
 */
BOOST_AUTO_TEST_CASE(mode_change_closes_popup)
{
    MutationFixture fixture;
    auto& sk = fixture.sk;
    sk.addFromString("PEPTIDE2{G}$$$$V2.0");
    auto* mol_model = sk.m_mol_model;
    mol_model->selectAll();

    auto* popup = open_monomer_db_popup(sk);
    BOOST_TEST(popup->isMutationMode());
    // a change that keeps a selection doesn't close the popup
    mol_model->select({mol_model->getMol()->getAtomWithIdx(0)}, {}, {}, {}, {},
                      SelectMode::DESELECT);
    BOOST_TEST(sk.m_sketcher_model->hasActiveSelection());
    BOOST_TEST(popup->isVisible());
    // clearing the selection does
    mol_model->clearSelection();
    BOOST_TEST(!popup->isVisible());

    popup = open_monomer_db_popup(sk);
    BOOST_TEST(!popup->isMutationMode());
    mol_model->selectAll();
    BOOST_TEST(!popup->isVisible());
}

/**
 * While the popup is open, Ctrl+A selects the search text rather than the
 * canvas contents, and does nothing in the table
 */
BOOST_AUTO_TEST_CASE(select_all_in_popup)
{
    MutationFixture fixture;
    auto& sk = fixture.sk;
    sk.m_mol_model->clearSelection();
    sk.show();
    sk.activateWindow();
    BOOST_REQUIRE(QTest::qWaitForWindowActive(&sk));

    auto* popup = open_monomer_db_popup(sk);
    auto* search_box = get_search_box(*popup);
    QTest::keyClicks(search_box, "pep");
    search_box->deselect();
    QTest::keyClick(search_box, Qt::Key_A, Qt::ControlModifier);
    BOOST_TEST(search_box->selectedText() == "pep");
    BOOST_TEST(!sk.m_sketcher_model->hasActiveSelection());
    BOOST_TEST(popup->isVisible());

    auto* table = get_monomer_db_table(*popup);
    table->setFocus();
    QTest::keyClick(table, Qt::Key_A, Qt::ControlModifier);
    BOOST_TEST(get_selected_rows(*table).empty());
    BOOST_TEST(!sk.m_sketcher_model->hasActiveSelection());
    BOOST_TEST(popup->isVisible());
    popup->hide();
}

/**
 * The popup is placed below or above its anchor, shrinking its table to fit
 * the available area, and clamped to it. If even its minimum layout doesn't
 * fit, its top-left corner stays inside the area.
 */
BOOST_AUTO_TEST_CASE(bounds)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    const auto full_size = popup->sizeHint();
    const auto minimum_size = popup->minimumSizeHint();
    BOOST_REQUIRE(minimum_size.height() < full_size.height());
    const QSize anchor_size(40, 30);
    const QRect bounds(0, 0, 1000, 800);

    // below the anchor, at full size
    QRect anchor({100, 100}, anchor_size);
    popup->updateGeometryToFit(anchor, bounds);
    BOOST_TEST(popup->size() == full_size);
    BOOST_TEST(popup->pos() == QPoint(100, anchor.bottom() + 1));

    // above the anchor, at full size
    anchor.moveTop(bounds.bottom() - 50);
    popup->updateGeometryToFit(anchor, bounds);
    BOOST_TEST(popup->size() == full_size);
    BOOST_TEST(popup->pos() == QPoint(100, anchor.top() - full_size.height()));

    // clamped at the right edge
    anchor.moveTopLeft({bounds.right() - 20, 100});
    popup->updateGeometryToFit(anchor, bounds);
    BOOST_TEST(popup->x() + popup->width() - 1 == bounds.right());

    // neither side has room, so the table shrinks to fit the side with more
    // room
    anchor = QRect({100, 50}, anchor_size);
    // bounds with the given room below the anchor, and less room above it
    auto bounds_with_room_below = [&anchor](int room_below) {
        return QRect(0, 0, 1000, anchor.bottom() + room_below + 1);
    };
    const int room_below = full_size.height() - 40;
    BOOST_REQUIRE(room_below > minimum_size.height());
    BOOST_REQUIRE(room_below > anchor.top());
    const int table_height = table->height();
    popup->updateGeometryToFit(anchor, bounds_with_room_below(room_below));
    BOOST_TEST(popup->height() == room_below);
    BOOST_TEST(popup->width() == full_size.width());
    BOOST_TEST(popup->pos() == QPoint(100, anchor.bottom() + 1));
    QApplication::processEvents();
    BOOST_TEST(table->height() < table_height);
    // the other controls stay accessible
    for (const auto* name :
         {"monomer_db_close_btn", "monomer_db_search_le",
          "monomer_db_close_link", "monomer_db_match_count_lbl"}) {
        auto* child = popup->findChild<QWidget*>(name);
        BOOST_TEST_CONTEXT(name)
        {
            BOOST_TEST(popup->rect().contains(
                QRect(child->mapTo(popup, QPoint(0, 0)), child->size())));
        }
    }

    // the available area shrinks while the popup is open, down to the
    // minimum layout
    const int smaller_room_below = minimum_size.height() + 10;
    const auto smaller_bounds = bounds_with_room_below(smaller_room_below);
    popup->updateGeometryToFit(anchor, smaller_bounds);
    BOOST_TEST(popup->height() == smaller_room_below);
    BOOST_TEST(popup->geometry().bottom() == smaller_bounds.bottom());

    // an area smaller than the minimum layout keeps the top-left corner inside
    // it, with the overflow to the right and bottom
    const QRect tiny_bounds(50, 60, minimum_size.width() - 100,
                            minimum_size.height() - 100);
    anchor = QRect({tiny_bounds.right() - 10, tiny_bounds.bottom() - 10},
                   anchor_size);
    popup->updateGeometryToFit(anchor, tiny_bounds);
    BOOST_TEST(popup->size() == minimum_size);
    BOOST_TEST(popup->pos() == tiny_bounds.topLeft());
    popup->hide();
}

/**
 * Click the filter button and return the filter popup, which must then be
 * visible
 */
static MonomerDatabaseFilterPopup*
open_filter_popup(MonomerDatabasePopup& popup)
{
    find_popup_child<QToolButton>(popup, "monomer_db_filter_btn")->click();
    auto* filter_popup = find_popup_child<MonomerDatabaseFilterPopup>(
        popup, "monomer_db_filter_popup");
    BOOST_REQUIRE(filter_popup->isVisible());
    return filter_popup;
}

/**
 * Check or uncheck a polymer type in the filter popup as the user would
 */
static void set_polymer_type_checked(MonomerDatabaseFilterPopup& filter_popup,
                                     const QString& text, bool checked)
{
    auto* combo =
        filter_popup.findChild<MultiSelectComboBox*>("polymer_type_combo");
    BOOST_REQUIRE(combo != nullptr);
    const int row = combo->findText(text);
    BOOST_REQUIRE(row >= 0);
    combo->setItemChecked(row, checked);
}

/**
 * Choose an item of one of the filter popup's single-select combo boxes as the
 * user would
 */
static void choose_filter_item(MonomerDatabaseFilterPopup& filter_popup,
                               const char* combo_name, const QString& text)
{
    auto* combo = filter_popup.findChild<QComboBox*>(combo_name);
    BOOST_REQUIRE(combo != nullptr);
    const int row = combo->findText(text);
    BOOST_REQUIRE_MESSAGE(row >= 0, text.toStdString());
    combo->setCurrentIndex(row);
}

/**
 * @return Whether the filter button shows the "filtering" (orange) icon
 */
static bool is_filter_button_filtering(MonomerDatabasePopup& popup)
{
    auto* button =
        find_popup_child<QToolButton>(popup, "monomer_db_filter_btn");
    return button->styleSheet() ==
           MONOMER_DATABASE_FILTER_BUTTON_FILTERING_STYLE;
}

/**
 * The filter popup opens from the filter button with each page's default
 * filter, and Clear All produces the most permissive filter on either page
 * while keeping the search text and sort. The filter button is orange
 * whenever any filter is applied.
 */
BOOST_AUTO_TEST_CASE(filter_defaults_and_clear_all)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    struct Case {
        bool amino_acid;
        QList<QVariant> checked_polymer_types;
        bool monomer_type_enabled;
    };
    for (const auto& [amino_acid, checked_polymer_types, monomer_type_enabled] :
         {Case{true,
               {static_cast<unsigned>(MonomerDbPolymerType::PEPTIDE),
                static_cast<unsigned>(MonomerDbPolymerType::CHEM)},
               true},
          Case{false,
               {static_cast<unsigned>(MonomerDbPolymerType::NUCLEIC_ACID),
                static_cast<unsigned>(MonomerDbPolymerType::CHEM)},
               false}}) {
        BOOST_TEST_CONTEXT("amino_acid=" << amino_acid)
        {
            fixture.clickPageButton(amino_acid);
            auto* popup = fixture.openPopup();
            BOOST_TEST(is_filter_button_filtering(*popup));
            QTest::keyClicks(get_search_box(*popup), "e");
            auto* table = get_monomer_db_table(*popup);
            table->sortByColumn(static_cast<int>(MonomerDbColumn::NAME),
                                Qt::DescendingOrder);
            const auto sort = popup->getSort();

            auto* filter_popup = open_filter_popup(*popup);
            auto* polymer_type_combo =
                filter_popup->findChild<MultiSelectComboBox*>(
                    "polymer_type_combo");
            BOOST_CHECK(polymer_type_combo->checkedData() ==
                        checked_polymer_types);
            BOOST_TEST(filter_popup->findChild<QComboBox*>("monomer_type_combo")
                           ->isEnabled() == monomer_type_enabled);
            BOOST_TEST(
                filter_popup->findChild<QComboBox*>("natural_analog_combo")
                    ->isEnabled());
            auto* filters_applied_label =
                filter_popup->findChild<QLabel*>("filters_applied_lbl");
            BOOST_TEST(filters_applied_label->text() == "1 filter applied");

            find_popup_child<QToolButton>(*popup, "filter_clear_all_link")
                ->click();
            BOOST_CHECK(popup->getFilter() == MonomerDatabaseFilter());
            BOOST_TEST(!is_filter_button_filtering(*popup));
            BOOST_TEST(polymer_type_combo->getDisplayText() == "Any");
            BOOST_TEST(!polymer_type_combo->font().bold());
            BOOST_TEST(filters_applied_label->text() == "No filters applied");
            BOOST_TEST(get_search_box(*popup)->text() == "e");
            BOOST_CHECK(popup->getSort() == sort);
            // every polymer type is shown, so the Type column is too
            BOOST_TEST(!table->isColumnHidden(
                static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));
            filter_popup->hide();
            popup->hide();
        }
    }
}

/**
 * Every filter change re-filters the table live. With {Peptide, Structure},
 * Terminal or a specific analog hides the structure rows and None keeps them.
 * The Type column hides when the filter allows a single polymer type, and no
 * polymer types empties the table.
 */
BOOST_AUTO_TEST_CASE(live_filtering)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    auto* filter_popup = open_filter_popup(*popup);
    auto check_symbols = [&](const std::vector<std::string>& symbols) {
        BOOST_TEST(get_monomer_db_symbols(*table) == symbols,
                   boost::test_tools::per_element());
        BOOST_TEST(get_match_count_text(*popup).toStdString() ==
                   std::to_string(symbols.size()) +
                       (symbols.size() == 1 ? " match" : " matches"));
    };
    check_symbols({"ac", "am", "chem1", "pepA", "pepX", "Zeta"});

    choose_filter_item(*filter_popup, "monomer_type_combo", "Terminal");
    check_symbols({"ac", "am"});
    choose_filter_item(*filter_popup, "monomer_type_combo", "Backbone");
    check_symbols({"pepA", "pepX", "Zeta"});
    choose_filter_item(*filter_popup, "monomer_type_combo", "Any");
    choose_filter_item(*filter_popup, "natural_analog_combo", "A");
    check_symbols({"pepA"});
    choose_filter_item(*filter_popup, "natural_analog_combo", "Any");
    check_symbols({"pepA", "Zeta"});
    choose_filter_item(*filter_popup, "natural_analog_combo", "None");
    check_symbols({"ac", "am", "chem1", "pepX"});
    BOOST_TEST(!table->isColumnHidden(
        static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));

    // a single polymer type hides the Type column
    set_polymer_type_checked(*filter_popup, "Structure", false);
    check_symbols({"ac", "am", "pepX"});
    BOOST_TEST(
        table->isColumnHidden(static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));
    BOOST_TEST(is_filter_button_filtering(*popup));

    // no polymer types shows nothing, and rechecking Peptide restores the
    // remembered analog
    set_polymer_type_checked(*filter_popup, "Peptide", false);
    check_symbols({});
    BOOST_TEST(is_filter_button_filtering(*popup));
    set_polymer_type_checked(*filter_popup, "Peptide", true);
    check_symbols({"ac", "am", "pepX"});

    // every polymer type shows everything and makes the button normal
    set_polymer_type_checked(*filter_popup, "Structure", true);
    set_polymer_type_checked(*filter_popup, "Nucleic acid", true);
    check_symbols({"ac", "am", "chem1", "pepA", "pepX", "rnaA", "Zeta"});
    BOOST_TEST(!table->isColumnHidden(
        static_cast<int>(MonomerDbColumn::POLYMER_TYPE)));
    BOOST_TEST(!is_filter_button_filtering(*popup));
    filter_popup->hide();
    popup->hide();
}

/**
 * Values hidden by applicability survive closing and reopening the filter
 * popup, and a page-default reset forgets them
 */
BOOST_AUTO_TEST_CASE(remembered_values_across_reopening)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto* popup = fixture.openPopup();
    auto* table = get_monomer_db_table(*popup);
    auto* filter_popup = open_filter_popup(*popup);
    set_polymer_type_checked(*filter_popup, "Structure", false);
    choose_filter_item(*filter_popup, "monomer_type_combo", "Backbone");
    choose_filter_item(*filter_popup, "natural_analog_combo", "A");
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"pepA"}),
               boost::test_tools::per_element());
    set_polymer_type_checked(*filter_popup, "Nucleic acid", true);
    BOOST_TEST(get_monomer_db_symbols(*table).size() == 6);
    filter_popup->hide();

    // the reopened popup shows the disabled combos' effective values
    filter_popup = open_filter_popup(*popup);
    auto* monomer_type_combo =
        filter_popup->findChild<QComboBox*>("monomer_type_combo");
    auto* analog_combo =
        filter_popup->findChild<QComboBox*>("natural_analog_combo");
    BOOST_TEST(monomer_type_combo->currentText() == "Any");
    BOOST_TEST(!monomer_type_combo->isEnabled());
    BOOST_TEST(analog_combo->currentText() == "Any or None");
    BOOST_TEST(!analog_combo->isEnabled());
    // restoring Peptide only restores Backbone and A
    set_polymer_type_checked(*filter_popup, "Nucleic acid", false);
    BOOST_TEST(monomer_type_combo->currentText() == "Backbone");
    BOOST_TEST(monomer_type_combo->font().bold());
    BOOST_TEST(analog_combo->currentText() == "A");
    BOOST_TEST(analog_combo->font().bold());
    BOOST_TEST(get_monomer_db_symbols(*table) ==
                   (std::vector<std::string>{"pepA"}),
               boost::test_tools::per_element());
    set_polymer_type_checked(*filter_popup, "Nucleic acid", true);
    filter_popup->hide();
    popup->hide();

    // a page-default reset forgets them. (Opening the popup on the other page
    // makes that the page at close, so that the popup resets when it opens on
    // this page again.)
    fixture.clickPageButton(false);
    fixture.openPopup()->hide();
    fixture.clickPageButton(true);
    popup = fixture.openPopup();
    BOOST_CHECK(popup->getFilter() ==
                MonomerDatabaseFilter(MonomerDbPolymerType::PEPTIDE |
                                      MonomerDbPolymerType::CHEM));
    popup->hide();
}

/**
 * Hiding the Monomer Database popup, e.g. with the load link, also hides the
 * filter popup, and reopening on the same page keeps the filter
 */
BOOST_AUTO_TEST_CASE(load_link_keeps_filter)
{
    PopupFixture fixture;
    auto* popup = fixture.openPopup();
    auto* filter_popup = open_filter_popup(*popup);
    choose_filter_item(*filter_popup, "monomer_type_combo", "Terminal");
    const auto filter = popup->getFilter();
    BOOST_REQUIRE(filter.getMonomerType() == MonomerTypeFilter::TERMINAL);

    get_load_link(*popup)->click();
    BOOST_TEST(!popup->isVisible());
    BOOST_TEST(!filter_popup->isVisible());
    get_file_chooser(fixture.widget)->reject();
    QApplication::processEvents();

    popup = fixture.openPopup();
    BOOST_CHECK(popup->getFilter() == filter);
    BOOST_TEST(get_monomer_db_symbols(*get_monomer_db_table(*popup)) ==
                   (std::vector<std::string>{"ac", "am"}),
               boost::test_tools::per_element());
    popup->hide();
}

/**
 * After enabling Nucleic acid on the amino acid page, choosing an RNA row
 * switches to the nucleic acid page while the popup stays open and anchored.
 * Closing and reopening it there keeps the filter and search text.
 */
BOOST_AUTO_TEST_CASE(cross_page_row_choice)
{
    PopupFixture fixture;
    load_custom_monomers(TEST_MONOMERS_JSON, NUM_TEST_MONOMERS);
    auto& model = fixture.model;
    auto* popup = fixture.openPopup();
    const auto popup_pos = popup->pos();
    QTest::keyClicks(get_search_box(*popup), "a");
    auto* filter_popup = open_filter_popup(*popup);
    set_polymer_type_checked(*filter_popup, "Nucleic acid", true);
    filter_popup->hide();
    const auto filter = popup->getFilter();

    click_monomer_db_row(*get_monomer_db_table(*popup), "rnaA");
    BOOST_CHECK(model.getMonomerToolType() == MonomerToolType::NUCLEIC_ACID);
    BOOST_CHECK(model.getMonomerDbMonomer() ==
                (MonomerID{"rnaA", ChainType::RNA}));
    BOOST_TEST(popup->isVisible());
    BOOST_TEST(popup->pos() == popup_pos);
    popup->hide();

    popup = fixture.openPopup();
    BOOST_CHECK(popup->getFilter() == filter);
    BOOST_TEST(get_search_box(*popup)->text() == "a");
    popup->hide();
}

/**
 * @return The widget after the given one that Tab moves the focus to, i.e. the
 * next visible, enabled widget in the focus chain that accepts Tab focus
 */
static QWidget* get_next_tab_stop(QWidget* widget)
{
    for (auto* next = widget->nextInFocusChain(); next != widget;
         next = next->nextInFocusChain()) {
        if ((next->focusPolicy() & Qt::TabFocus) && next->isVisible() &&
            next->isEnabled()) {
            return next;
        }
    }
    return nullptr;
}

/**
 * Tab moves the focus from the search box to the filter button, then to the
 * table and then to the links
 */
BOOST_AUTO_TEST_CASE(tab_order)
{
    PopupFixture fixture;
    auto* popup = fixture.openPopup();
    QWidget* widget = get_search_box(*popup);
    for (const auto* name : {"monomer_db_filter_btn", "monomer_db_table",
                             "monomer_db_load_link", "monomer_db_close_link"}) {
        BOOST_TEST_CONTEXT(name)
        {
            widget = get_next_tab_stop(widget);
            BOOST_REQUIRE(widget != nullptr);
            BOOST_TEST(widget->objectName().toStdString() == name);
        }
    }
    popup->hide();
}

/**
 * With a canvas selection under Select or Move/Rotate, navigating to an RNA
 * row doesn't switch pages, while committing it does, keeping the draw tool
 * and each page's last monomeric tool
 */
BOOST_AUTO_TEST_CASE(cross_page_mutation)
{
    for (const auto draw_tool : {DrawTool::SELECT, DrawTool::MOVE_ROTATE}) {
        BOOST_TEST_CONTEXT("draw_tool=" << static_cast<int>(draw_tool))
        {
            MutationFixture fixture;
            auto& sk = fixture.sk;
            auto* sketcher_model = sk.m_sketcher_model;
            sketcher_model->setValue(ModelKey::DRAW_TOOL, draw_tool);
            BOOST_REQUIRE(sketcher_model->hasActiveSelection());
            const auto amino_acid_tool = sketcher_model->getAminoAcidTool();
            const auto aa_restore = sketcher_model->getDrawToolToRestore(
                MonomerToolType::AMINO_ACID);
            const auto na_restore = sketcher_model->getDrawToolToRestore(
                MonomerToolType::NUCLEIC_ACID);

            auto* popup = open_monomer_db_popup(sk);
            BOOST_REQUIRE(popup->isMutationMode());
            auto* filter_popup = open_filter_popup(*popup);
            set_polymer_type_checked(*filter_popup, "Nucleic acid", true);
            filter_popup->hide();
            auto* table = get_monomer_db_table(*popup);
            table->setFocus();
            while (get_current_symbol(*table) != "rnaA") {
                const auto previous = get_current_symbol(*table);
                QTest::keyClick(table, Qt::Key_Down);
                BOOST_REQUIRE(get_current_symbol(*table) != previous);
            }
            BOOST_CHECK(sketcher_model->getMonomerToolType() ==
                        MonomerToolType::AMINO_ACID);

            QTest::keyClick(table, Qt::Key_Return);
            BOOST_CHECK(sketcher_model->getMonomerToolType() ==
                        MonomerToolType::NUCLEIC_ACID);
            BOOST_CHECK(sketcher_model->getDrawTool() == draw_tool);
            BOOST_CHECK(sketcher_model->getAminoAcidTool() == amino_acid_tool);
            BOOST_CHECK(sketcher_model->getDrawToolToRestore(
                            MonomerToolType::AMINO_ACID) == aa_restore);
            BOOST_CHECK(sketcher_model->getDrawToolToRestore(
                            MonomerToolType::NUCLEIC_ACID) == na_restore);
            // the selected peptide can't become an RNA monomer
            BOOST_TEST(fixture.getHelm() == "PEPTIDE1{A}$$$$V2.0");
            BOOST_TEST(popup->isVisible());
            popup->hide();
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
