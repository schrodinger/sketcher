#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include <boost/signals2/connection.hpp>
#include <QMetaObject>
#include <QPointer>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/monomer_database_filter.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"

class QModelIndex;
class QRect;

namespace Ui
{
class MonomerDatabasePopup;
}

namespace schrodinger
{
namespace sketcher
{

enum class MonomerToolType;
class MonomerDatabaseFilterPopup;
class MonomerDatabaseTableModel;
class SketcherModel;

/**
 * A popup for searching the monomer database and choosing a monomer from it.
 *
 * The table lists the core monomers with no natural analog plus the custom
 * monomers that don't shadow a core monomer with a natural analog (see
 * get_monomer_database_table_monomers()), filtered and sorted without
 * re-reading the database.
 *
 * The popup has two modes, which are fixed while it's open:
 *  - In tool mode (no canvas selection), changing rows with the mouse or
 *    keyboard switches the scene tool to the row's monomer.
 *  - In mutation mode (canvas selection), clicking on a row mutates all
 *    selected monomers to match the clicked monomer.  With the keyboard, the up
 *    and down arrows move the table's current row but don't trigger a mutation.
 *    Hitting Enter, Return, or Space mutate all selected monomers to match the
 *    monomer of the table's current row.
 *
 * The popup closes itself if a canvas selection appears or disappears while
 * it's open.
 *
 * A table row is selected only when the model's active draw tool is the
 * current page's database monomer, there's no canvas selection, and that
 * monomer is in the table.
 */
class SKETCHER_API MonomerDatabasePopup : public SketcherViewWithWasmOutlineFix
{
    Q_OBJECT
  public:
    /**
     * @param model The model to read the tool state, side bar page and canvas
     * selection from. It must outlive the popup.
     * @param anchor The widget that the popup is always shown next to, i.e. the
     * side bar's monomer database button, however the popup is opened
     * @param parent The parent widget, which is also used as the parent of the
     * file chooser and dialogs opened by the load database link
     */
    MonomerDatabasePopup(const SketcherModel* model, QWidget* anchor,
                         QWidget* parent = nullptr);
    ~MonomerDatabasePopup();

    /**
     * Show the popup next to the anchor given at construction. If this is the
     * first opening, or the model's side bar page differs from the page when
     * the popup last closed, the search text is cleared and the filter is
     * reset to the page's default. Otherwise, the search text, filter and sort
     * are kept.
     *
     * Opening the popup never chooses a monomer.
     */
    void showPopup();

    /**
     * Select the table row of the model's active database monomer, if there
     * is one (see the class documentation). Otherwise clear the table
     * selection, but keep the current index. This never chooses a monomer.
     */
    void syncSelectionToModel();

    /**
     * Set whether the Load Database... or Update Database... link is visible
     */
    void setLoadDatabaseVisible(bool visible);

    /**
     * @return The current filter
     */
    const MonomerDatabaseFilter& getFilter() const;

    /**
     * @return The current sort
     */
    const MonomerDbSort& getSort() const;

    /**
     * @return Whether the popup was opened in mutation mode, i.e. with a canvas
     * selection
     */
    bool isMutationMode() const;

    /**
     * Resize the popup to fit next to the anchor inside the bounds, and move it
     * there. The table shrinks to keep the popup inside the bounds, down to
     * the popup's minimum size. If the minimum size doesn't fit, the popup's
     * top-left corner is kept inside the bounds.
     *
     * The popup calls this when it opens and whenever its screen's available
     * geometry changes.
     *
     * @param anchor_rect The anchor's global geometry
     * @param bounds The area to fit the popup inside, in global coordinates
     */
    void updateGeometryToFit(const QRect& anchor_rect, const QRect& bounds);

  signals:
    /**
     * Emitted when the user chooses a monomer from the table
     */
    void monomerChosen(const rdkit_extensions::MonomerID& monomer);

  protected:
    /**
     * Overridden QWidget method: record the page at close, however the popup
     * closes
     */
    void hideEvent(QHideEvent* event) override;

  private:
    /**
     * Re-read the table's monomers from the database, then update the table
     */
    void refreshDatabase();

    /**
     * Re-run the filter and sort, then update the table and the widgets that
     * depend on the filter (the filter button and the Type column)
     */
    void refilter();

    /**
     * Respond to the user changing the filter in the filter popup
     */
    void onFilterChanged(const MonomerDatabaseFilter& filter);

    /**
     * Open the filter popup, initialized from the current filter
     */
    void onFilterButtonClicked();

    /**
     * Apply the search text to the filtered rows and reset the table model,
     * preserving the current row (by monomer ID) when it's still in the table,
     * then synchronize the selection with the model
     */
    void updateTableRows();

    /**
     * Respond to the user clicking a column header
     */
    void onSortRequested(const MonomerDbSort& sort);

    /**
     * Set the header's sort indicator to match the stored sort without
     * requesting another sort
     */
    void updateSortIndicator();

    /**
     * In tool mode, choose the monomer of a row that the user selected
     */
    void onTableSelectionChanged();

    /**
     * In mutation mode, choose the monomer of the clicked row
     */
    void onTableClicked(const QModelIndex& index);

    /**
     * Choose the monomer of the given row
     */
    void chooseRow(int row);

    /**
     * Close the popup if a canvas selection appeared or disappeared
     */
    void onModelSelectionChanged();

    /**
     * Close the popup and ask the user for a monomer database to load
     */
    void onLoadLinkClicked();

    /**
     * Refit the popup next to its anchor after the screen's available geometry
     * changes
     */
    void onAvailableGeometryChanged();

    /**
     * @return The table row of the model's active database monomer, if the
     * table should select it
     */
    std::optional<int> getModelSelectedRow() const;

    std::unique_ptr<Ui::MonomerDatabasePopup> m_ui;
    const SketcherModel* m_model;
    MonomerDatabaseTableModel* m_table_model = nullptr;
    MonomerDatabaseFilterPopup* m_filter_popup = nullptr;
    // The table's source list, which is shared with the table model
    std::shared_ptr<const std::vector<rdkit_extensions::MonomerInfo>>
        m_all_monomers;
    // Indices into m_all_monomers that pass the filter, in sorted order
    std::vector<std::size_t> m_filtered_rows;
    MonomerDatabaseFilter m_filter;
    MonomerDbSort m_sort;
    // The side bar page when the popup last closed, or nullopt if it never
    // has
    std::optional<MonomerToolType> m_page_at_close;
    bool m_mutation_mode = false;
    // whether the table is being updated programmatically, during which
    // selection changes must not choose a monomer
    bool m_updating_table = false;
    // The side bar's monomer database button
    QPointer<QWidget> m_anchor;
    std::vector<QMetaObject::Connection> m_screen_connections;
    boost::signals2::scoped_connection m_database_connection;
};

} // namespace sketcher
} // namespace schrodinger
