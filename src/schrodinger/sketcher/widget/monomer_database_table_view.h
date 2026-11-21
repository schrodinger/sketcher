#pragma once

#include <QItemSelectionModel>
#include <QTableView>

#include "schrodinger/sketcher/definitions.h"

class QKeyEvent;
class QPaintEvent;

namespace schrodinger
{
namespace sketcher
{

/**
 * The table of the Monomer Database popup. It
 *  - selects single rows, except in mutation mode, where the user can move the
 *    current index but never select a row
 *  - highlights the current row when it isn't the selected row
 *  - reports Enter, Return and Space as requests to commit the current row
 *  - shows a message when it has no rows
 *  - has fixed row heights and no vertical header
 */
class SKETCHER_API MonomerDatabaseTableView : public QTableView
{
    Q_OBJECT
  public:
    MonomerDatabaseTableView(QWidget* parent = nullptr);

    /**
     * Set whether the table is in mutation mode, in which mouse and keyboard
     * interaction only moves the current index and never selects a row.
     * Programmatic selection is unaffected.
     */
    void setMutationMode(bool mutation_mode);

    /**
     * @return Whether the table is in mutation mode
     */
    bool isMutationMode() const;

    /**
     * Overridden QAbstractScrollArea method: room for the header and
     * PREFERRED_VISIBLE_ROWS rows
     */
    QSize sizeHint() const override;

    /**
     * Overridden QAbstractScrollArea method: room for the header and
     * MIN_VISIBLE_ROWS rows
     */
    QSize minimumSizeHint() const override;

    /// The number of rows that the table shows when there's enough room
    static constexpr int PREFERRED_VISIBLE_ROWS = 10;
    /// The number of rows that the table shows when it has been shrunk to fit
    /// in a small area
    static constexpr int MIN_VISIBLE_ROWS = 2;

  signals:
    /**
     * Emitted when the user presses Enter, Return or Space while the table has
     * focus and a valid current index. Auto-repeated key presses are ignored.
     * @param row The current row
     */
    void commitRequested(int row);

  protected:
    QItemSelectionModel::SelectionFlags
    selectionCommand(const QModelIndex& index,
                     const QEvent* event = nullptr) const override;
    void keyPressEvent(QKeyEvent* event) override;
    void currentChanged(const QModelIndex& current,
                        const QModelIndex& previous) override;
    void paintEvent(QPaintEvent* event) override;

  private:
    /**
     * @return The height of the table with room for the given number of rows
     */
    int calcHeightForRows(int num_rows) const;

    /**
     * Repaint the whole of the given row, if it's valid
     */
    void updateRow(int row);

    bool m_mutation_mode = false;
};

} // namespace sketcher
} // namespace schrodinger
