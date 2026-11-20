#define BOOST_TEST_MODULE test_monomer_database_table_model

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <QSignalSpy>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/widget/monomer_database_table_model.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

using schrodinger::rdkit_extensions::ChainType;
using schrodinger::rdkit_extensions::MonomerID;
using schrodinger::rdkit_extensions::MonomerInfo;

namespace schrodinger
{
namespace sketcher
{

MAKE_ENUM_LOGGABLE(MonomerDbColumn)

static MonomerInfo make_monomer(const std::string& symbol,
                                const std::string& polymer_type,
                                std::optional<std::string> natural_analog,
                                std::optional<std::string> name)
{
    MonomerInfo info;
    info.symbol = symbol;
    info.polymer_type = polymer_type;
    info.natural_analog = std::move(natural_analog);
    info.name = std::move(name);
    return info;
}

static std::shared_ptr<const std::vector<MonomerInfo>> get_monomers()
{
    return std::make_shared<const std::vector<MonomerInfo>>(
        std::vector<MonomerInfo>{
            make_monomer("dC", "PEPTIDE", "C", "D-Cysteine"),
            make_monomer("ac", "PEPTIDE", "X", "Acetyl"),
            make_monomer("Eth", "CHEM", std::nullopt, std::nullopt),
            make_monomer("5mC", "RNA", "C", "5-methylcytosine"),
            make_monomer("C", "RNA", "N", "Not really cytosine"),
        });
}

static QString get_display_text(const MonomerDatabaseTableModel& model, int row,
                                MonomerDbColumn column)
{
    return model.data(model.index(row, static_cast<int>(column))).toString();
}

BOOST_AUTO_TEST_CASE(test_empty_model)
{
    MonomerDatabaseTableModel model;
    BOOST_TEST(model.rowCount() == 0);
    BOOST_TEST(model.columnCount() == 4);
}

BOOST_AUTO_TEST_CASE(test_set_monomers_and_display_data)
{
    MonomerDatabaseTableModel model;
    QSignalSpy about_to_reset_spy(&model,
                                  &QAbstractItemModel::modelAboutToBeReset);
    QSignalSpy reset_spy(&model, &QAbstractItemModel::modelReset);
    auto monomers = get_monomers();
    model.setMonomers(monomers, {3, 0, 2, 1});
    BOOST_TEST(about_to_reset_spy.count() == 1);
    BOOST_TEST(reset_spy.count() == 1);

    BOOST_TEST(model.rowCount() == 4);
    BOOST_TEST(model.columnCount() == 4);
    // nothing is nested under a row
    BOOST_TEST(model.rowCount(model.index(0, 0)) == 0);
    BOOST_TEST(model.columnCount(model.index(0, 0)) == 0);

    std::vector<std::vector<QString>> expected{
        {"5mC", "5-methylcytosine", "Nucleic acid", "C"},
        {"dC", "D-Cysteine", "Peptide", "C"},
        // a CHEM monomer with no name or analog
        {"Eth", "", "Structure", ""},
        // a None analog (X) is blank
        {"ac", "Acetyl", "Peptide", ""},
    };
    for (int row = 0; row < model.rowCount(); ++row) {
        for (int col = 0; col < model.columnCount(); ++col) {
            BOOST_TEST_CONTEXT("row " << row << ", column " << col)
            {
                BOOST_TEST(get_display_text(
                               model, row, static_cast<MonomerDbColumn>(col)) ==
                           expected[row][col]);
            }
        }
    }

    // the analog column is centered, and the others use the default alignment
    auto analog_index =
        model.index(0, static_cast<int>(MonomerDbColumn::NATURAL_ANALOG));
    BOOST_TEST((model.data(analog_index, Qt::TextAlignmentRole)
                    .value<Qt::Alignment>() == Qt::AlignCenter));
    auto symbol_index =
        model.index(0, static_cast<int>(MonomerDbColumn::SYMBOL));
    BOOST_TEST(!model.data(symbol_index, Qt::TextAlignmentRole).isValid());
    BOOST_TEST(!model.data(symbol_index, Qt::ToolTipRole).isValid());
    BOOST_TEST(!model.data(QModelIndex()).isValid());

    // the model shares the list rather than copying it
    BOOST_TEST(&model.getMonomerInfo(1) == &(*monomers)[0]);

    // a None analog for a nucleic acid (N) is blank too
    model.setMonomers(monomers, {4});
    BOOST_TEST(reset_spy.count() == 2);
    BOOST_TEST(model.rowCount() == 1);
    BOOST_TEST(get_display_text(model, 0, MonomerDbColumn::NATURAL_ANALOG) ==
               "");

    model.setMonomers(monomers, {});
    BOOST_TEST(reset_spy.count() == 3);
    BOOST_TEST(model.rowCount() == 0);
}

BOOST_AUTO_TEST_CASE(test_header_labels)
{
    MonomerDatabaseTableModel model;
    std::vector<QString> expected{"Symbol", "Name", "Type", "Analog"};
    for (int col = 0; col < model.columnCount(); ++col) {
        BOOST_TEST(model.headerData(col, Qt::Horizontal).toString() ==
                   expected[col]);
    }
    BOOST_TEST(!model.headerData(4, Qt::Horizontal).isValid());
    BOOST_TEST(!model.headerData(0, Qt::Horizontal, Qt::ToolTipRole).isValid());
}

BOOST_AUTO_TEST_CASE(test_row_and_id_lookups)
{
    MonomerDatabaseTableModel model;
    model.setMonomers(get_monomers(), {3, 0, 2});

    BOOST_TEST((model.getMonomerID(0) == MonomerID{"5mC", ChainType::RNA}));
    BOOST_TEST((model.getMonomerID(1) == MonomerID{"dC", ChainType::PEPTIDE}));
    BOOST_TEST((model.getMonomerID(2) == MonomerID{"Eth", ChainType::CHEM}));
    BOOST_TEST(*model.getMonomerInfo(2).symbol == "Eth");

    for (int row = 0; row < model.rowCount(); ++row) {
        BOOST_TEST(model.getRow(model.getMonomerID(row)).value_or(-1) == row);
    }
    // the same symbol with a different polymer type
    BOOST_TEST(!model.getRow({"dC", ChainType::RNA}).has_value());
    // a monomer that's in the list but isn't displayed
    BOOST_TEST(!model.getRow({"ac", ChainType::PEPTIDE}).has_value());
    BOOST_TEST(!model.getRow({"missing", ChainType::PEPTIDE}).has_value());
}

BOOST_AUTO_TEST_CASE(test_sort_signal)
{
    MonomerDatabaseTableModel model;
    model.setMonomers(get_monomers(), {0, 1, 2});
    QSignalSpy sort_spy(&model, &MonomerDatabaseTableModel::sortRequested);
    QSignalSpy reset_spy(&model, &QAbstractItemModel::modelReset);
    QSignalSpy layout_spy(&model, &QAbstractItemModel::layoutChanged);

    model.sort(static_cast<int>(MonomerDbColumn::POLYMER_TYPE),
               Qt::DescendingOrder);
    BOOST_REQUIRE(sort_spy.count() == 1);
    auto sort = sort_spy.at(0).at(0).value<MonomerDbSort>();
    BOOST_TEST(sort.column == MonomerDbColumn::POLYMER_TYPE);
    BOOST_TEST(!sort.ascending);

    model.sort(static_cast<int>(MonomerDbColumn::NAME), Qt::AscendingOrder);
    BOOST_REQUIRE(sort_spy.count() == 2);
    sort = sort_spy.at(1).at(0).value<MonomerDbSort>();
    BOOST_TEST((sort == MonomerDbSort{MonomerDbColumn::NAME, true}));

    // the model doesn't sort itself
    BOOST_TEST(reset_spy.count() == 0);
    BOOST_TEST(layout_spy.count() == 0);
    BOOST_TEST(get_display_text(model, 0, MonomerDbColumn::SYMBOL) == "dC");

    // invalid columns (e.g. clearing the sort indicator) are ignored
    model.sort(-1, Qt::AscendingOrder);
    model.sort(4, Qt::AscendingOrder);
    BOOST_TEST(sort_spy.count() == 2);
}

} // namespace sketcher
} // namespace schrodinger
