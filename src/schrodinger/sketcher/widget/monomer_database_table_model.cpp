#include "schrodinger/sketcher/widget/monomer_database_table_model.h"

#include <utility>

#include <QString>

using schrodinger::rdkit_extensions::MonomerID;
using schrodinger::rdkit_extensions::MonomerInfo;

namespace schrodinger
{
namespace sketcher
{

MonomerDatabaseTableModel::MonomerDatabaseTableModel(QObject* parent) :
    QAbstractTableModel(parent)
{
}

void MonomerDatabaseTableModel::setMonomers(
    std::shared_ptr<const std::vector<MonomerInfo>> all_monomers,
    std::vector<std::size_t> rows)
{
    beginResetModel();
    m_all_monomers = std::move(all_monomers);
    m_rows = std::move(rows);
    endResetModel();
}

int MonomerDatabaseTableModel::rowCount(const QModelIndex& parent) const
{
    if (parent.isValid()) {
        return 0;
    }
    return static_cast<int>(m_rows.size());
}

int MonomerDatabaseTableModel::columnCount(const QModelIndex& parent) const
{
    if (parent.isValid()) {
        return 0;
    }
    return static_cast<int>(MONOMER_DB_COLUMN_LABELS.size());
}

QVariant MonomerDatabaseTableModel::data(const QModelIndex& index,
                                         int role) const
{
    if (!index.isValid() || index.row() >= rowCount() ||
        index.column() >= columnCount()) {
        return QVariant();
    }
    auto column = static_cast<MonomerDbColumn>(index.column());
    if (role == Qt::TextAlignmentRole &&
        column == MonomerDbColumn::NATURAL_ANALOG) {
        return QVariant::fromValue(Qt::Alignment(Qt::AlignCenter));
    } else if (role != Qt::DisplayRole) {
        return QVariant();
    }

    const auto& info = getMonomerInfo(index.row());
    switch (column) {
        case MonomerDbColumn::SYMBOL:
            return QString::fromStdString(info.symbol.value_or(""));
        case MonomerDbColumn::NAME:
            return QString::fromStdString(info.name.value_or(""));
        case MonomerDbColumn::POLYMER_TYPE:
            return QString::fromStdString(get_polymer_type_display_text(
                get_monomer_db_polymer_type(info)));
        case MonomerDbColumn::NATURAL_ANALOG:
            return QString::fromStdString(
                get_natural_analog_display_text(info));
    }
    return QVariant();
}

QVariant MonomerDatabaseTableModel::headerData(int section,
                                               Qt::Orientation orientation,
                                               int role) const
{
    if (orientation != Qt::Horizontal) {
        return QAbstractTableModel::headerData(section, orientation, role);
    }
    if (role == Qt::DisplayRole && section >= 0 && section < columnCount()) {
        auto label = MONOMER_DB_COLUMN_LABELS[section].second;
        return QString::fromUtf8(label.data(),
                                 static_cast<qsizetype>(label.size()));
    }
    return QVariant();
}

void MonomerDatabaseTableModel::sort(int column, Qt::SortOrder order)
{
    if (column < 0 || column >= columnCount()) {
        return;
    }
    emit sortRequested(
        {static_cast<MonomerDbColumn>(column), order == Qt::AscendingOrder});
}

MonomerID MonomerDatabaseTableModel::getMonomerID(int row) const
{
    // The table only contains filtered rows, which always have a valid polymer
    // type, so this never throws
    return get_monomer_id(getMonomerInfo(row)).value();
}

std::optional<int>
MonomerDatabaseTableModel::getRow(const MonomerID& monomer_id) const
{
    for (int row = 0; row < rowCount(); ++row) {
        // check the symbol first to avoid parsing every row's polymer type
        if (getMonomerInfo(row).symbol == monomer_id.symbol &&
            getMonomerID(row) == monomer_id) {
            return row;
        }
    }
    return std::nullopt;
}

const MonomerInfo& MonomerDatabaseTableModel::getMonomerInfo(int row) const
{
    return m_all_monomers->at(m_rows.at(row));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/monomer_database_table_model.moc"
