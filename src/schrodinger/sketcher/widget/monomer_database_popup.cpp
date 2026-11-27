#include "schrodinger/sketcher/widget/monomer_database_popup.h"

#include <algorithm>
#include <string>
#include <utility>

#include <QHeaderView>
#include <QHideEvent>
#include <QIcon>
#include <QItemSelectionModel>
#include <QLineEdit>
#include <QMetaObject>
#include <QScreen>
#include <QSignalBlocker>

#include "schrodinger/sketcher/dialog/custom_title_bar.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_monomer_database_popup.h"
#include "schrodinger/sketcher/widget/monomer_database_table_model.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

using schrodinger::rdkit_extensions::MonomerDatabase;
using schrodinger::rdkit_extensions::MonomerInfo;

namespace schrodinger
{
namespace sketcher
{

static const QString SEARCH_ICON_PATH = ":icons/search.svg";
static const QString LOAD_DATABASE_TEXT = "Load Database...";
static const QString UPDATE_DATABASE_TEXT = "Update Database...";

/**
 * @return The filter to use when the popup opens on the given page
 */
static MonomerDatabaseFilter get_page_default_filter(MonomerToolType page)
{
    return MonomerDatabaseFilter((page == MonomerToolType::AMINO_ACID
                                      ? MonomerDbPolymerType::PEPTIDE
                                      : MonomerDbPolymerType::NUCLEIC_ACID) |
                                 MonomerDbPolymerType::CHEM);
}

/**
 * @return The text for the match count label
 */
static QString get_match_count_text(int num_matches)
{
    return QString("%1 %2")
        .arg(num_matches)
        .arg(num_matches == 1 ? "match" : "matches");
}

MonomerDatabasePopup::MonomerDatabasePopup(const SketcherModel* model,
                                           QWidget* anchor, QWidget* parent) :
    SketcherViewWithWasmOutlineFix(parent),
    m_model(model),
    m_anchor(anchor)
{
    setObjectName("monomer_db_popup");
    m_ui.reset(new Ui::MonomerDatabasePopup());
    m_ui->setupUi(this);
    setWindowFlags(Qt::Popup);
    setStyleSheet(MONOMER_DATABASE_POPUP_STYLE);
    // A popup never has a native title bar
    auto* title_bar =
        new CustomTitleBar(windowTitle(), this, TitleBarType::POPUP);
    title_bar->setObjectName("monomer_db_title_bar");
    title_bar->m_title_label->setObjectName("monomer_db_title_lbl");
    auto* close_button = title_bar->getCloseButton();
    close_button->setObjectName("monomer_db_close_btn");
    m_ui->main_layout->insertWidget(0, title_bar);
    m_ui->monomer_db_load_link->setStyleSheet(TEXT_LINK_STYLE);
    m_ui->monomer_db_close_link->setStyleSheet(TEXT_LINK_STYLE);
    m_ui->monomer_db_search_le->addAction(QIcon(SEARCH_ICON_PATH),
                                          QLineEdit::LeadingPosition);

    m_table_model = new MonomerDatabaseTableModel(this);
    auto* table = m_ui->monomer_db_table;
    table->setModel(m_table_model);
    auto* header = table->horizontalHeader();
    header->setSectionResizeMode(QHeaderView::Interactive);
    header->setSectionResizeMode(static_cast<int>(MonomerDbColumn::NAME),
                                 QHeaderView::Stretch);
    table->setColumnWidth(static_cast<int>(MonomerDbColumn::SYMBOL), 90);
    table->setColumnWidth(static_cast<int>(MonomerDbColumn::POLYMER_TYPE), 85);
    table->setColumnWidth(static_cast<int>(MonomerDbColumn::NATURAL_ANALOG),
                          55);
    // Enabling sorting immediately sorts by the header's sort indicator, so
    // set the indicator to the stored sort first
    updateSortIndicator();
    table->setSortingEnabled(true);

    connect(m_table_model, &MonomerDatabaseTableModel::sortRequested, this,
            &MonomerDatabasePopup::onSortRequested);
    connect(table->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MonomerDatabasePopup::onTableSelectionChanged);
    connect(table, &MonomerDatabaseTableView::clicked, this,
            &MonomerDatabasePopup::onTableClicked);
    connect(table, &MonomerDatabaseTableView::commitRequested, this,
            &MonomerDatabasePopup::chooseRow);
    connect(m_ui->monomer_db_search_le, &QLineEdit::textChanged, this,
            &MonomerDatabasePopup::updateTableRows);
    connect(close_button, &QToolButton::clicked, this,
            &MonomerDatabasePopup::hide);
    connect(m_ui->monomer_db_close_link, &QToolButton::clicked, this,
            &MonomerDatabasePopup::hide);
    connect(m_ui->monomer_db_load_link, &QToolButton::clicked, this,
            &MonomerDatabasePopup::onLoadLinkClicked);
    connect(m_model, &SketcherModel::selectionChanged, this,
            &MonomerDatabasePopup::onModelSelectionChanged);

    // Read the database directly rather than relying on any other database
    // listener (such as SketcherModel) having run first
    m_database_connection =
        MonomerDatabase::instance().subscribeToChanges([this]() {
            QMetaObject::invokeMethod(
                this, [this]() { refreshDatabase(); }, Qt::AutoConnection);
        });
    m_filter = get_page_default_filter(m_model->getMonomerToolType());
    refreshDatabase();
}

MonomerDatabasePopup::~MonomerDatabasePopup() = default;

void MonomerDatabasePopup::refreshDatabase()
{
    const auto& db = MonomerDatabase::instance();
    auto custom_monomers = db.getCustomMonomers();
    const bool has_custom_monomers = !custom_monomers.empty();
    m_all_monomers = std::make_shared<const std::vector<MonomerInfo>>(
        get_monomer_database_table_monomers(std::move(custom_monomers),
                                            db.getCoreMonomers()));
    m_ui->monomer_db_load_link->setText(
        has_custom_monomers ? UPDATE_DATABASE_TEXT : LOAD_DATABASE_TEXT);
    refilter();
}

void MonomerDatabasePopup::refilter()
{
    m_filtered_rows =
        filter_and_sort_monomers(*m_all_monomers, m_filter, m_sort);
    // The Type column is only useful when more than one polymer type is shown
    m_ui->monomer_db_table->setColumnHidden(
        static_cast<int>(MonomerDbColumn::POLYMER_TYPE),
        count_flags(m_filter.getPolymerTypes()) == 1);
    updateTableRows();
}

void MonomerDatabasePopup::updateTableRows()
{
    auto* table = m_ui->monomer_db_table;
    auto* selection_model = table->selectionModel();
    // Remember the current row by ID, since its row number may change
    std::optional<rdkit_extensions::MonomerID> current_id;
    if (auto current = table->currentIndex(); current.isValid()) {
        current_id = m_table_model->getMonomerID(current.row());
    }
    m_updating_table = true;
    auto rows =
        filter_by_search_text(*m_all_monomers, m_filtered_rows,
                              m_ui->monomer_db_search_le->text().toStdString());
    m_table_model->setMonomers(m_all_monomers, std::move(rows));
    if (current_id.has_value()) {
        if (auto row = m_table_model->getRow(*current_id)) {
            // Restore the browsing row. If the monomer disappeared or was
            // filtered out, the current index stays cleared rather than moving
            // to another monomer.
            QSignalBlocker blocker(selection_model);
            selection_model->setCurrentIndex(m_table_model->index(*row, 0),
                                             QItemSelectionModel::NoUpdate);
        }
    }
    m_updating_table = false;
    m_ui->monomer_db_match_count_lbl->setText(
        get_match_count_text(m_table_model->rowCount()));
    // A row selected by the model takes precedence over the browsing row
    syncSelectionToModel();
}

std::optional<int> MonomerDatabasePopup::getModelSelectedRow() const
{
    if (m_model->getDrawTool() != DrawTool::MONOMER_DB_MONOMER ||
        m_model->hasActiveSelection()) {
        return std::nullopt;
    }
    // A monomer that isn't in the table (e.g. because the database changed and
    // the model hasn't been updated yet) selects nothing
    return m_table_model->getRow(m_model->getMonomerDbMonomer());
}

void MonomerDatabasePopup::syncSelectionToModel()
{
    auto* table = m_ui->monomer_db_table;
    auto* selection_model = table->selectionModel();
    {
        // Block the signals so that synchronizing never chooses a monomer
        QSignalBlocker blocker(selection_model);
        if (auto row = getModelSelectedRow()) {
            auto index = m_table_model->index(*row, 0);
            selection_model->setCurrentIndex(
                index, QItemSelectionModel::ClearAndSelect |
                           QItemSelectionModel::Rows);
            table->scrollTo(index);
        } else {
            // Unlike clear(), this keeps the current index
            selection_model->clearSelection();
        }
    }
    // The view didn't see the blocked signals, so repaint it explicitly
    table->viewport()->update();
}

void MonomerDatabasePopup::showPopup()
{
    if (!m_anchor) {
        return;
    }
    const auto page = m_model->getMonomerToolType();
    if (!m_page_at_close.has_value() || *m_page_at_close != page) {
        {
            QSignalBlocker blocker(m_ui->monomer_db_search_le);
            m_ui->monomer_db_search_le->clear();
        }
        m_filter = get_page_default_filter(page);
    }
    m_mutation_mode = m_model->hasActiveSelection();
    m_ui->monomer_db_table->setMutationMode(m_mutation_mode);
    // Don't keep highlighting the row that was current when the popup last
    // closed. If the model has a row selected, refiltering makes it current
    // again.
    m_ui->monomer_db_table->selectionModel()->clearCurrentIndex();
    refilter();

    // On WASM builds, each canvas is its own screen, so this also notices the
    // canvas resizing
    auto* screen = m_anchor->screen();
    for (auto signal :
         {&QScreen::availableGeometryChanged, &QScreen::geometryChanged}) {
        m_screen_connections.push_back(
            connect(screen, signal, this,
                    &MonomerDatabasePopup::onAvailableGeometryChanged));
    }
    updateGeometryToFit(get_global_geometry(m_anchor),
                        get_popup_bounds(m_anchor));
    // Showing a popup gives the focus to the widget in it that last had the
    // focus. Make that the search box before showing, since a table that gets
    // the focus without a current index makes its first row current.
    m_ui->monomer_db_search_le->setFocus(Qt::PopupFocusReason);
    show_popup_next_to_widget(this, m_anchor);
}

void MonomerDatabasePopup::updateGeometryToFit(const QRect& anchor_rect,
                                               const QRect& bounds)
{
    ensurePolished();
    layout()->activate();
    const auto full_size = sizeHint();
    const auto minimum_size = minimumSizeHint();
    // Shrink the popup (and so the table, which is the only part that can
    // shrink) to fit above or below the anchor, but no smaller than the
    // minimum layout
    const int height =
        std::max(std::min(full_size.height(),
                          calc_max_popup_height(anchor_rect, bounds)),
                 minimum_size.height());
    const int width = std::max(std::min(full_size.width(), bounds.width()),
                               minimum_size.width());
    resize(width, height);
    move(calc_popup_position(size(), anchor_rect, bounds));
}

void MonomerDatabasePopup::onAvailableGeometryChanged()
{
    if (isVisible() && m_anchor) {
        updateGeometryToFit(get_global_geometry(m_anchor),
                            get_popup_bounds(m_anchor));
    }
}

void MonomerDatabasePopup::hideEvent(QHideEvent* event)
{
    SketcherViewWithWasmOutlineFix::hideEvent(event);
    m_page_at_close = m_model->getMonomerToolType();
    for (const auto& connection : m_screen_connections) {
        disconnect(connection);
    }
    m_screen_connections.clear();
}

void MonomerDatabasePopup::onSortRequested(const MonomerDbSort& sort)
{
    if (sort == m_sort) {
        return;
    }
    m_sort = sort;
    updateSortIndicator();
    refilter();
}

void MonomerDatabasePopup::updateSortIndicator()
{
    auto* header = m_ui->monomer_db_table->horizontalHeader();
    QSignalBlocker blocker(header);
    header->setSortIndicator(static_cast<int>(m_sort.column),
                             m_sort.ascending ? Qt::AscendingOrder
                                              : Qt::DescendingOrder);
}

void MonomerDatabasePopup::onTableSelectionChanged()
{
    if (m_updating_table || m_mutation_mode) {
        return;
    }
    const auto selected =
        m_ui->monomer_db_table->selectionModel()->selectedRows();
    if (selected.size() == 1) {
        chooseRow(selected.front().row());
    }
}

void MonomerDatabasePopup::onTableClicked(const QModelIndex& index)
{
    if (m_mutation_mode && index.isValid()) {
        chooseRow(index.row());
    }
}

void MonomerDatabasePopup::chooseRow(int row)
{
    if (row < 0 || row >= m_table_model->rowCount()) {
        return;
    }
    emit monomerChosen(m_table_model->getMonomerID(row));
}

void MonomerDatabasePopup::onModelSelectionChanged()
{
    // The mode is fixed while the popup is open. The scene may re-emit
    // selectionChanged without the selection changing (e.g. while rebuilding
    // items after a mutation), so only close if the mode would change.
    if (isVisible() && m_model->hasActiveSelection() != m_mutation_mode) {
        hide();
    }
}

void MonomerDatabasePopup::onLoadLinkClicked()
{
    // Close the popup before the file chooser opens. It stays closed
    // afterwards, however the load turns out. The chooser and result dialogs
    // are parented to the popup's parent, since the popup is hidden.
    hide();
    load_monomer_database_from_file(parentWidget());
}

void MonomerDatabasePopup::setLoadDatabaseVisible(bool visible)
{
    m_ui->monomer_db_load_link->setVisible(visible);
}

const MonomerDatabaseFilter& MonomerDatabasePopup::getFilter() const
{
    return m_filter;
}

const MonomerDbSort& MonomerDatabasePopup::getSort() const
{
    return m_sort;
}

bool MonomerDatabasePopup::isMutationMode() const
{
    return m_mutation_mode;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/monomer_database_popup.moc"
