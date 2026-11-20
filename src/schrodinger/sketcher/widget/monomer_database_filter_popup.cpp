#include "schrodinger/sketcher/widget/monomer_database_filter_popup.h"

#include <string>
#include <utility>
#include <vector>

#include <QAbstractItemView>
#include <QComboBox>
#include <QHideEvent>
#include <QList>
#include <QSignalBlocker>
#include <QString>
#include <QVariant>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_monomer_database_filter_popup.h"
#include "schrodinger/sketcher/widget/multi_select_combo_box.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

static const std::vector<MonomerDbPolymerType> POLYMER_TYPE_ITEMS = {
    MonomerDbPolymerType::PEPTIDE, MonomerDbPolymerType::NUCLEIC_ACID,
    MonomerDbPolymerType::CHEM};

static const std::vector<std::pair<MonomerTypeFilter, QString>>
    MONOMER_TYPE_ITEMS = {{MonomerTypeFilter::ANY, "Any"},
                          {MonomerTypeFilter::BACKBONE, "Backbone"},
                          {MonomerTypeFilter::TERMINAL, "Terminal"}};

// The natural analog combo box's items that aren't a specific natural analog.
// The specific natural analogs follow them.
static const std::vector<std::pair<NaturalAnalogFilter, QString>>
    NATURAL_ANALOG_ITEMS = {{NaturalAnalogFilter::ANY_OR_NONE, "Any or None"},
                            {NaturalAnalogFilter::ANY, "Any"},
                            {NaturalAnalogFilter::NONE, "None"}};

/**
 * @return The text for the filter count label
 */
static QString get_filters_applied_text(int num_filters)
{
    if (num_filters == 0) {
        return "No filters applied";
    }
    return QString("%1 %2 applied")
        .arg(num_filters)
        .arg(num_filters == 1 ? "filter" : "filters");
}

/**
 * Make the combo box's text bold or not, without making its dropdown's items
 * bold. This must be called again after any items are added to the combo box.
 */
static void set_combo_bold(QComboBox* combo, bool bold)
{
    auto font = combo->font();
    font.setBold(false);
    // The dropdown draws an item using the item's Qt::FontRole data if there
    // is any, and the combo box's font otherwise
    for (int row = 0; row < combo->count(); ++row) {
        combo->setItemData(row, font, Qt::FontRole);
    }
    font.setBold(bold);
    combo->setFont(font);
}

MonomerDatabaseFilterPopup::MonomerDatabaseFilterPopup(QWidget* anchor,
                                                       QWidget* parent) :
    SketcherViewWithWasmOutlineFix(parent),
    m_anchor(anchor)
{
    setObjectName("monomer_db_filter_popup");
    m_ui.reset(new Ui::MonomerDatabaseFilterPopup());
    m_ui->setupUi(this);
    setWindowFlags(Qt::Popup);
    setStyleSheet(MONOMER_DATABASE_FILTER_POPUP_STYLE);
    m_ui->filter_clear_all_link->setStyleSheet(TEXT_LINK_STYLE);

    auto* polymer_type_combo = m_ui->polymer_type_combo;
    for (auto polymer_type : POLYMER_TYPE_ITEMS) {
        polymer_type_combo->addItem(
            QString::fromStdString(get_polymer_type_display_text(polymer_type)),
            static_cast<unsigned>(polymer_type));
    }
    for (const auto& [monomer_type, text] : MONOMER_TYPE_ITEMS) {
        m_ui->monomer_type_combo->addItem(text, static_cast<int>(monomer_type));
    }
    updateWidgets();

    connect(polymer_type_combo, &MultiSelectComboBox::checkedItemsChanged, this,
            &MonomerDatabaseFilterPopup::onPolymerTypesChanged);
    connect(m_ui->monomer_type_combo, &QComboBox::currentIndexChanged, this,
            &MonomerDatabaseFilterPopup::onMonomerTypeChanged);
    connect(m_ui->natural_analog_combo, &QComboBox::currentIndexChanged, this,
            &MonomerDatabaseFilterPopup::onNaturalAnalogChanged);
    connect(m_ui->filter_clear_all_link, &QToolButton::clicked, this,
            &MonomerDatabaseFilterPopup::onClearAllClicked);
    connect(m_ui->filter_close_btn, &QToolButton::clicked, this,
            &MonomerDatabaseFilterPopup::hide);
}

MonomerDatabaseFilterPopup::~MonomerDatabaseFilterPopup() = default;

void MonomerDatabaseFilterPopup::showPopup(const MonomerDatabaseFilter& filter)
{
    if (!m_anchor) {
        return;
    }
    m_filter = filter;
    updateWidgets();
    show_popup_next_to_widget(this, m_anchor);
}

void MonomerDatabaseFilterPopup::reposition()
{
    if (isVisible() && m_anchor) {
        move_popup_next_to_widget(this, m_anchor);
    }
}

const MonomerDatabaseFilter& MonomerDatabaseFilterPopup::getFilter() const
{
    return m_filter;
}

void MonomerDatabaseFilterPopup::hideEvent(QHideEvent* event)
{
    SketcherViewWithWasmOutlineFix::hideEvent(event);
    for (auto* combo : findChildren<QComboBox*>()) {
        combo->hidePopup();
    }
}

void MonomerDatabaseFilterPopup::onPolymerTypesChanged()
{
    auto polymer_types = MonomerDbPolymerType::NONE;
    for (const auto& data : m_ui->polymer_type_combo->checkedData()) {
        polymer_types =
            polymer_types | static_cast<MonomerDbPolymerType>(data.toUInt());
    }
    m_filter.setPolymerTypes(polymer_types);
    onFilterEdited();
}

void MonomerDatabaseFilterPopup::onMonomerTypeChanged()
{
    m_filter.setMonomerType(static_cast<MonomerTypeFilter>(
        m_ui->monomer_type_combo->currentData().toInt()));
    onFilterEdited();
}

void MonomerDatabaseFilterPopup::onNaturalAnalogChanged()
{
    if (!m_filter.isNaturalAnalogApplicable()) {
        // the combo box is disabled, so the user can't have changed it
        return;
    }
    auto* combo = m_ui->natural_analog_combo;
    const auto analog =
        static_cast<NaturalAnalogFilter>(combo->currentData().toInt());
    m_filter.setNaturalAnalog(analog, analog == NaturalAnalogFilter::SPECIFIC
                                          ? combo->currentText().toStdString()
                                          : std::string());
    onFilterEdited();
}

void MonomerDatabaseFilterPopup::onClearAllClicked()
{
    m_filter = MonomerDatabaseFilter();
    onFilterEdited();
}

void MonomerDatabaseFilterPopup::onFilterEdited()
{
    updateWidgets();
    emit filterChanged(m_filter);
}

void MonomerDatabaseFilterPopup::updateWidgets()
{
    {
        auto* combo = m_ui->polymer_type_combo;
        QSignalBlocker blocker(combo);
        QList<QVariant> checked_data;
        for (auto polymer_type : POLYMER_TYPE_ITEMS) {
            if (has_flag(m_filter.getPolymerTypes(), polymer_type)) {
                checked_data.append(static_cast<unsigned>(polymer_type));
            }
        }
        combo->setCheckedData(checked_data);
        set_combo_bold(combo, m_filter.isPolymerTypeFiltered());
        // the display text doesn't depend on the current item, so make sure
        // that the combo box repaints
        combo->update();
    }
    {
        auto* combo = m_ui->monomer_type_combo;
        QSignalBlocker blocker(combo);
        set_combo_box_data(combo, static_cast<int>(m_filter.getMonomerType()));
        combo->setEnabled(m_filter.isMonomerTypeApplicable());
        set_combo_bold(combo, m_filter.isMonomerTypeFiltered());
    }
    updateNaturalAnalogCombo();
    m_ui->filters_applied_lbl->setText(
        get_filters_applied_text(m_filter.getNumFiltersApplied()));
}

void MonomerDatabaseFilterPopup::updateNaturalAnalogCombo()
{
    auto* combo = m_ui->natural_analog_combo;
    QSignalBlocker blocker(combo);
    combo->clear();
    for (const auto& [analog, text] : NATURAL_ANALOG_ITEMS) {
        combo->addItem(text, static_cast<int>(analog));
    }
    for (const auto& choice : m_filter.getNaturalAnalogChoices()) {
        combo->addItem(QString::fromStdString(choice),
                       static_cast<int>(NaturalAnalogFilter::SPECIFIC));
    }
    const auto analog = m_filter.getNaturalAnalog();
    if (analog == NaturalAnalogFilter::SPECIFIC) {
        combo->setCurrentIndex(combo->findText(
            QString::fromStdString(m_filter.getSpecificNaturalAnalog())));
    } else {
        set_combo_box_data(combo, static_cast<int>(analog));
    }
    combo->setEnabled(m_filter.isNaturalAnalogApplicable());
    set_combo_bold(combo, m_filter.isNaturalAnalogFiltered());
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/monomer_database_filter_popup.moc"
