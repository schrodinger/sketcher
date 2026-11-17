#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include <QAbstractTableModel>
#include <QMetaType>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/monomer_database_filter.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A table model for the Monomer Database popup. The model is a thin adapter
 * over a shared list of monomers and the filtered and sorted row indices into
 * that list (see filter_and_sort_monomers() and filter_by_search_text()). It
 * contains no filtering or sorting logic of its own. Its columns are given by
 * MonomerDbColumn.
 */
class SKETCHER_API MonomerDatabaseTableModel : public QAbstractTableModel
{
    Q_OBJECT
  public:
    explicit MonomerDatabaseTableModel(QObject* parent = nullptr);

    /**
     * Replace the contents of the model.
     * @param all_monomers All monomers that `rows` may index into. The model
     * shares ownership of this list.
     * @param rows The indices into `all_monomers` of the monomers to display,
     * in display order
     */
    void setMonomers(
        std::shared_ptr<const std::vector<rdkit_extensions::MonomerInfo>>
            all_monomers,
        std::vector<std::size_t> rows);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    int columnCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index,
                  int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const override;

    /**
     * Request a sort by emitting sortRequested(). This method does not sort
     * the model itself; the owner is expected to re-sort the rows and call
     * setMonomers().
     */
    void sort(int column, Qt::SortOrder order = Qt::AscendingOrder) override;

    /**
     * @return the ID of the monomer displayed in the given row
     * @pre `row` is a valid row and its monomer has a valid polymer type
     */
    rdkit_extensions::MonomerID getMonomerID(int row) const;

    /**
     * @return the row that displays the given monomer, or std::nullopt if the
     * monomer isn't displayed
     */
    std::optional<int>
    getRow(const rdkit_extensions::MonomerID& monomer_id) const;

    /**
     * @return the monomer displayed in the given row
     * @pre `row` is a valid row
     */
    const rdkit_extensions::MonomerInfo& getMonomerInfo(int row) const;

  signals:
    /**
     * Emitted when the view requests a new sort order
     */
    void sortRequested(const MonomerDbSort& sort);

  private:
    std::shared_ptr<const std::vector<rdkit_extensions::MonomerInfo>>
        m_all_monomers;
    std::vector<std::size_t> m_rows;
};

} // namespace sketcher
} // namespace schrodinger

Q_DECLARE_METATYPE(schrodinger::sketcher::MonomerDbSort);
