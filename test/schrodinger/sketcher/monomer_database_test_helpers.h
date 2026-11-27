/**
 * @file
 * @brief Helpers for C++ tests that load a custom monomer database or use the
 *        Monomer Database popup.
 */

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <QAbstractItemModel>
#include <QApplication>
#include <QTableView>
#include <QTest>
#include <QToolButton>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/widget/monomer_database_popup.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Undo any changes to the monomer database when the test that owns this
 * object finishes, even if the test fails
 */
struct ResetMonomerDatabase {
    ~ResetMonomerDatabase()
    {
        rdkit_extensions::MonomerDatabase::instance().resetMonomerDefinitions();
    }
};

/**
 * Load the given JSON as the custom monomer database, replacing any existing
 * custom monomers, and require that every monomer loads without errors
 * @param json The monomer definitions
 * @param num_monomers The number of monomers in json
 */
inline void load_custom_monomers(std::string_view json,
                                 std::size_t num_monomers)
{
    auto result =
        rdkit_extensions::MonomerDatabase::instance().loadMonomersFromJson(
            json);
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), num_monomers);
}

/**
 * Click the Monomer Database button in the given widget (e.g. a
 * MonomerToolWidget or a SketcherWidget) and return the popup, which must
 * then be visible
 */
inline MonomerDatabasePopup* open_monomer_db_popup(QWidget& widget)
{
    auto* button = widget.findChild<QToolButton*>("monomer_db_btn");
    BOOST_REQUIRE(button != nullptr);
    button->click();
    auto* popup = widget.findChild<MonomerDatabasePopup*>("monomer_db_popup");
    BOOST_REQUIRE(popup != nullptr);
    BOOST_REQUIRE(popup->isVisible());
    return popup;
}

/**
 * @return The Monomer Database popup's table
 */
inline QTableView* get_monomer_db_table(MonomerDatabasePopup& popup)
{
    auto* table = popup.findChild<QTableView*>("monomer_db_table");
    BOOST_REQUIRE(table != nullptr);
    return table;
}

/**
 * @return The symbols in the table, in display order
 */
inline std::vector<std::string> get_monomer_db_symbols(const QTableView& table)
{
    std::vector<std::string> symbols;
    const auto* model = table.model();
    for (int row = 0; row < model->rowCount(); ++row) {
        symbols.push_back(model->index(row, 0).data().toString().toStdString());
    }
    return symbols;
}

/**
 * @return The row that shows the given symbol and polymer type text (e.g.
 * "Peptide"), or the first row with the symbol if no type is given
 */
inline std::optional<int>
find_monomer_db_row(const QTableView& table, const std::string& symbol,
                    const std::string& polymer_type = "")
{
    const auto* model = table.model();
    for (int row = 0; row < model->rowCount(); ++row) {
        if (model->index(row, 0).data().toString().toStdString() == symbol &&
            (polymer_type.empty() ||
             model->index(row, 2).data().toString().toStdString() ==
                 polymer_type)) {
            return row;
        }
    }
    return std::nullopt;
}

/**
 * Click the given row of the table with the mouse, scrolling it into view first
 */
inline void click_monomer_db_row(QTableView& table, int row)
{
    auto index = table.model()->index(row, 0);
    table.scrollTo(index);
    QTest::mouseClick(table.viewport(), Qt::LeftButton, Qt::NoModifier,
                      table.visualRect(index).center());
    QApplication::processEvents();
}

/**
 * Click the row that shows the given symbol, which must be in the table
 */
inline void click_monomer_db_row(QTableView& table, const std::string& symbol,
                                 const std::string& polymer_type = "")
{
    auto row = find_monomer_db_row(table, symbol, polymer_type);
    BOOST_REQUIRE_MESSAGE(row.has_value(), symbol);
    click_monomer_db_row(table, *row);
}

} // namespace sketcher
} // namespace schrodinger
