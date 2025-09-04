#include "schrodinger/sketcher/dialog/edit_atom_properties.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/StereoGroup.h>

#include <QProxyStyle>
#include <QPushButton>
#include <QStringList>
#include <QTimer>

#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/ui/ui_edit_atom_properties.h"
#include "schrodinger/sketcher/widget/periodic_table_widget.h"
#include "schrodinger/sketcher/widget/blankable_spin_box.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

Q_DECLARE_METATYPE(RDKit::StereoGroupType);
Q_DECLARE_METATYPE(schrodinger::sketcher::QueryType);
Q_DECLARE_METATYPE(schrodinger::sketcher::AtomQuery);
Q_DECLARE_METATYPE(schrodinger::sketcher::QueryAromaticity);
Q_DECLARE_METATYPE(schrodinger::sketcher::QueryCount);

namespace schrodinger
{
namespace sketcher
{

// these constants define the display text and the Qt.UserRole data for all
// combo boxes in the dialog
const std::vector<std::pair<QString, RDKit::StereoGroupType>>
    STEREO_COMBO_DATA = {
        {"ABS", RDKit::StereoGroupType::STEREO_ABSOLUTE},
        {"AND", RDKit::StereoGroupType::STEREO_AND},
        {"OR", RDKit::StereoGroupType::STEREO_OR},
};

const std::vector<std::pair<QString, QueryType>> QUERY_TYPE_COMBO_DATA = {
    {"Allowed List", QueryType::ALLOWED_LIST},
    {"Not Allowed List", QueryType::NOT_ALLOWED_LIST},
    {"Wildcard", QueryType::WILDCARD},
    {"Specific Element", QueryType::SPECIFIC_ELEMENT},
    {"R-Group", QueryType::RGROUP},
    {"SMARTS", QueryType::SMARTS},
};

const std::vector<std::pair<QString, AtomQuery>> WILDCARD_COMBO_DATA = {
    {"Any heavy atom (A)", AtomQuery::A}, {"Heteroatom (Q)", AtomQuery::Q},
    {"Metal (M)", AtomQuery::M},          {"Halogen (X)", AtomQuery::X},
    {"Any or H (AH)", AtomQuery::AH},     {"Hetero or H (QH)", AtomQuery::QH},
    {"Metal or H (MH)", AtomQuery::MH},   {"Halogen or H (XH)", AtomQuery::XH},
};

const std::vector<std::pair<QString, QueryAromaticity>> AROMATICITY_COMBO_DATA =
    {
        {"(any)", QueryAromaticity::ANY},
        {"Aromatic (a)", QueryAromaticity::AROMATIC},
        {"Aliphatic (A)", QueryAromaticity::ALIPHATIC},
};

const std::vector<std::pair<QString, QueryCount>>
    ANY_POSITIVE_EXACTLY_COMBO_DATA = {
        {"(any)", QueryCount::ANY},
        {">0", QueryCount::POSITIVE},
        {"exactly", QueryCount::EXACTLY},
};

const std::vector<std::pair<QString, QueryCount>> ANY_EXACTLY_COMBO_DATA = {
    {"(any)", QueryCount::ANY},
    {"exactly", QueryCount::EXACTLY},
};

// Tooltips that apply to more than one widget. Other tooltips are set in the
// ui file
const QString ISOTOPE_TOOLTIP = "Isotope number";
const QString CHARGE_TOOLTIP = "Formal charge";
const QString UNPAIRED_TOOLTIP =
    "Number of unpaired electrons - defines the radical "
    "character of the system";
const QString STEREO_TOOLTIP =
    "Chiral center label for enhanced stereochemistry - ABS (absolute) is the "
    "default; AND / OR require a group identifier";

/**
 * Load the specified data into the combo box
 * @param combo The combo box to populate
 * @param data The data to load, which should be formatted as a vector of
 * (display text, Qt.UserRole data)
 */
template <typename T> static void
load_data_into_combo_box(QComboBox* const combo,
                         const std::vector<std::pair<QString, T>>& data)
{
    for (auto [text, val] : data) {
        combo->addItem(text, QVariant::fromValue(val));
    }
}

/**
 * Connect signals such that the spin box is only visible when
 * QueryCount::EXACTLY is selected in the combo box (as the Qt.UserRole data)
 */
static void
show_spin_box_when_exactly_is_selected_in_combo_box(QComboBox* const combo,
                                                    QSpinBox* const sb)
{
    combo->connect(combo, &QComboBox::currentIndexChanged, sb,
                   [sb, combo](const int index) {
                       auto variant = combo->itemData(index);
                       auto data = variant.value<QueryCount>();
                       sb->setVisible(data == QueryCount::EXACTLY);
                   });
}

/**
 * Connect signals such that the spin box is only visible when anything other
 * than EnhancedStereoType::ABS is selected in the combo box (as the Qt.UserRole
 * data)
 */
static void
hide_spin_box_when_abs_is_selected_in_combo_box(QComboBox* const combo,
                                                QSpinBox* const sb)
{
    combo->connect(combo, &QComboBox::currentIndexChanged, sb,
                   [sb, combo](const int index) {
                       auto variant = combo->itemData(index);
                       auto data = variant.value<RDKit::StereoGroupType>();
                       sb->setVisible(data !=
                                      RDKit::StereoGroupType::STEREO_ABSOLUTE);
                   });
}

EditAtomPropertiesDialog::EditAtomPropertiesDialog(const RDKit::Atom* atom,
                                                   MolModel* const mol_model,
                                                   QWidget* parent) :
    ModalDialog(parent),
    m_atom(atom),
    m_mol_model(mol_model)
{
    ui.reset(new Ui::EditAtomPropertiesDialog());
    setupDialogUI(*ui);

    m_update_buttons_enabled_timer = new QTimer(this);
    m_update_buttons_enabled_timer->setSingleShot(true);
    m_update_buttons_enabled_timer->setInterval(0);
    connect(m_update_buttons_enabled_timer, &QTimer::timeout, this,
            &EditAtomPropertiesDialog::updateButtonsEnabled);

    // connect combo boxes that update whether other widgets are visible
    show_spin_box_when_exactly_is_selected_in_combo_box(ui->total_h_combo,
                                                        ui->total_h_sb);
    show_spin_box_when_exactly_is_selected_in_combo_box(ui->ring_count_combo,
                                                        ui->ring_count_sb);
    show_spin_box_when_exactly_is_selected_in_combo_box(
        ui->ring_bond_count_combo, ui->ring_bond_count_sb);
    hide_spin_box_when_abs_is_selected_in_combo_box(ui->atom_stereo_combo,
                                                    ui->atom_stereo_sb);
    hide_spin_box_when_abs_is_selected_in_combo_box(ui->query_stereo_combo,
                                                    ui->query_stereo_sb);
    connect(ui->query_type_combo, &QComboBox::currentIndexChanged, this,
            &EditAtomPropertiesDialog::onQueryTypeComboBoxChanged);

    // populate the combo boxes
    load_data_into_combo_box(ui->atom_stereo_combo, STEREO_COMBO_DATA);
    load_data_into_combo_box(ui->query_stereo_combo, STEREO_COMBO_DATA);
    load_data_into_combo_box(ui->query_type_combo, QUERY_TYPE_COMBO_DATA);
    load_data_into_combo_box(ui->wildcard_combo, WILDCARD_COMBO_DATA);
    load_data_into_combo_box(ui->aromaticity_combo, AROMATICITY_COMBO_DATA);
    load_data_into_combo_box(ui->total_h_combo,
                             ANY_POSITIVE_EXACTLY_COMBO_DATA);
    load_data_into_combo_box(ui->ring_count_combo,
                             ANY_POSITIVE_EXACTLY_COMBO_DATA);
    load_data_into_combo_box(ui->ring_bond_count_combo, ANY_EXACTLY_COMBO_DATA);

    // set the tooltips that would have to be duplicated if we set them in the
    // ui file
    ui->atom_isotope_lbl->setToolTip(ISOTOPE_TOOLTIP);
    ui->query_isotope_lbl->setToolTip(ISOTOPE_TOOLTIP);
    ui->atom_charge_lbl->setToolTip(CHARGE_TOOLTIP);
    ui->query_charge_lbl->setToolTip(CHARGE_TOOLTIP);
    ui->atom_unpaired_lbl->setToolTip(UNPAIRED_TOOLTIP);
    ui->query_unpaired_lbl->setToolTip(UNPAIRED_TOOLTIP);
    ui->atom_stereo_lbl->setToolTip(STEREO_TOOLTIP);
    ui->query_stereo_lbl->setToolTip(STEREO_TOOLTIP);

    // set the spin box ranges
    ui->atom_isotope_sb->setRange(0, MAX_ISOTOPE_VALUE);
    ui->query_isotope_sb->setRange(0, MAX_ISOTOPE_VALUE);
    ui->atom_charge_sb->setRange(-ATOM_CHARGE_LIMIT, ATOM_CHARGE_LIMIT);
    ui->query_charge_sb->setRange(-ATOM_CHARGE_LIMIT, ATOM_CHARGE_LIMIT);
    ui->atom_unpaired_sb->setRange(MIN_UNPAIRED_E, MAX_UNPAIRED_E);
    ui->query_unpaired_sb->setRange(MIN_UNPAIRED_E, MAX_UNPAIRED_E);
    ui->atom_stereo_sb->setRange(1, MAX_STEREO_GROUP_ID);
    ui->query_stereo_sb->setRange(1, MAX_STEREO_GROUP_ID);
    ui->total_h_sb->setRange(0, MAX_QUERY_TOTAL_H);
    ui->num_connections_sb->setRange(0, MAX_QUERY_NUM_CONNECTIONS);
    ui->ring_count_sb->setRange(0, MAX_QUERY_RING_COUNT);
    ui->ring_bond_count_sb->setRange(0, MAX_QUERY_RING_BOND_COUNT);
    ui->smallest_ring_size_sb->setRange(0, MAX_QUERY_SMALLEST_RING_SIZE);

    // add validators to the line edits used for entering elements
    auto* element_validator =
        new ElementValidator(/*allow_list = */ false, this);
    auto* element_list_validator =
        new ElementValidator(/*allow_list = */ true, this);
    ui->atom_element_le->setValidator(element_validator);
    ui->query_element_le->setValidator(element_validator);
    ui->element_list_le->setValidator(element_list_validator);

    // hook up the atom/query radio buttons
    ui->set_as_group->setId(
        ui->set_as_atom_rb,
        ui->atom_query_stacked_wdg->indexOf(ui->edit_atom_page));
    ui->set_as_group->setId(
        ui->set_as_query_rb,
        ui->atom_query_stacked_wdg->indexOf(ui->edit_query_page));
    connect(ui->set_as_group, &QButtonGroup::idClicked,
            ui->atom_query_stacked_wdg, &QStackedWidget::setCurrentIndex);
    connect(ui->atom_query_stacked_wdg, &QStackedWidget::currentChanged, this,
            &EditAtomPropertiesDialog::onAtomOrQueryToggled);

    // update whether the OK and reset buttons are enabled any time anything
    // changes in the dialog. We use a timer to prevent updateButtonsEnabled
    // from running every time loadProperties updates each individual widget
    connect_input_widgets_to_timer(this, m_update_buttons_enabled_timer);

    // Set up the periodic table buttons and widgets
    std::vector<std::pair<ToolButtonWithPopup*,
                          void (EditAtomPropertiesDialog::*)(Element)>>
        periodic_table_mappings = {
            {ui->atom_periodic_table_btn,
             &EditAtomPropertiesDialog::onAtomPeriodicTableElementSelected},
            {ui->query_periodic_table_btn,
             &EditAtomPropertiesDialog::onQueryPeriodicTableElementSelected}};
    for (auto [btn, slot_func] : periodic_table_mappings) {
        auto* periodic_table_widget = new PeriodicTableWidget(this);
        periodic_table_widget->setWindowFlags(Qt::Popup);
        // we'll close the widget manually in the slots since we want it to stay
        // open when we're selecting into a list
        periodic_table_widget->setCloseOnClick(false);
        btn->setPopupDelay(0);
        btn->showPopupIndicator(false);
        btn->setPopupWidget(periodic_table_widget);
        // We still want to leave space for the query periodic table button even
        // when it's hidden. Otherwise, all of the other widgets will be moved
        // slightly whenever the button's visibility is toggled. (The atom
        // periodic table button is never hidden, so this property has no effect
        // on that button.)
        auto size_policy = btn->sizePolicy();
        size_policy.setRetainSizeWhenHidden(true);
        btn->setSizePolicy(size_policy);
        connect(periodic_table_widget, &PeriodicTableWidget::elementSelected,
                this, slot_func);
    }

    changeAllLineEditClearButtonIcons();

    auto* reset_btn = ui->buttonBox->button(QDialogButtonBox::Reset);
    connect(reset_btn, &QPushButton::clicked, this,
            &EditAtomPropertiesDialog::reset);
    // load in settings from the atom
    reset();
}

EditAtomPropertiesDialog::~EditAtomPropertiesDialog() = default;

void EditAtomPropertiesDialog::switchToAllowedList()
{
    ui->set_as_query_rb->click();
    set_combo_box_data(ui->query_type_combo, QueryType::ALLOWED_LIST);
}

void EditAtomPropertiesDialog::reset()
{
    auto props = read_properties_from_atom(m_atom);
    loadProperties(props);
    updateButtonsEnabled();
}

void EditAtomPropertiesDialog::changeAllLineEditClearButtonIcons()
{
    QIcon new_icon(LINE_EDIT_CLEAR_ICON_PATH);
    for (auto* line_edit : findChildren<QLineEdit*>()) {
        auto* clear_btn = line_edit->findChild<QToolButton*>();
        if (clear_btn != nullptr) {
            clear_btn->setIcon(new_icon);
        }
    }
}

/**
 * Convert the given text to an Element enum value
 * @throw InvalidAtomPropertyError if the text does not specify a valid element
 * (e.g "Qz", "banana", or the empty string)
 */
static Element get_element_from_text(const QString& qtext)
{
    if (qtext.isEmpty()) {
        throw InvalidAtomPropertyError("No element specified.");
    }
    auto text = qtext.toStdString();
    int element_num;
    try {
        // note that symbol_to_atomic_num is case insensitive
        element_num = symbol_to_atomic_number(text);
    } catch (std::runtime_error&) {
        throw InvalidAtomPropertyError(text + " is not a valid element.");
    }
    return Element(element_num);
}

/**
 * Split the given string that contains a list of element names into a list of
 * strings that contain element names (e.g. convert "C, N, O" into ["C", "N",
 * "O"])
 */
static QStringList split_element_list(const QString& text)
{
    return text.split(QRegularExpression(R"(\W+)"), Qt::SkipEmptyParts);
}

/**
 * Convert the given text into a vector of Element enum values
 * @throw InvalidAtomPropertyError if the text does not specify any elements
 * (e.g. the string is empty) or if the text specifies any invalid or duplicate
 * element (e.g. "C, N, Qz")
 */
static std::vector<Element> get_elements_from_list(const QString& text)
{
    auto split_element_text = split_element_list(text);
    std::vector<Element> elements;
    std::transform(split_element_text.begin(), split_element_text.end(),
                   std::back_inserter(elements), get_element_from_text);

    if (elements.empty()) {
        throw InvalidAtomPropertyError("No elements specified.");
    }
    //  throw when elements contains duplicates
    auto elements_copy = elements;
    std::sort(elements_copy.begin(), elements_copy.end());
    auto last = std::unique(elements_copy.begin(), elements_copy.end());
    if (last != elements_copy.end()) {
        throw InvalidAtomPropertyError("Duplicate elements specified.");
    }
    return elements;
}

/**
 * Set the enhanced stereo properties based on the values selected in the combo
 * box and spin box
 */
static void set_enhanced_stereo_properties(const QComboBox* const combo,
                                           const QSpinBox* const sb,
                                           AbstractAtomProperties& props)
{
    if (combo->isEnabled()) {
        EnhancedStereo enhanced_stereo;
        enhanced_stereo.setType(
            combo->currentData().value<RDKit::StereoGroupType>());
        enhanced_stereo.setGroupId(!sb->isHidden() ? sb->value() : 0);
        props.enhanced_stereo = enhanced_stereo;
    }
}

/**
 * Load the given enhanced stereo properties into the specified widgets
 */
static void set_enhanced_stereo_widgets(QLabel* const label,
                                        QComboBox* const combo,
                                        QSpinBox* const sb,
                                        const AbstractAtomProperties& props)
{
    bool enable = props.enhanced_stereo.has_value();
    // we want the combo box to show absolute stereo if it's disabled
    EnhancedStereo enhanced_stereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 1);
    if (enable) {
        enhanced_stereo = *props.enhanced_stereo;
    }
    label->setEnabled(enable);
    combo->setEnabled(enable);
    set_combo_box_data(combo, enhanced_stereo.type());
    // setting the combo box data will automatically show or hide the spin box
    sb->setValue(enhanced_stereo.groupId());
}

void EditAtomPropertiesDialog::accept()
{
    auto props = getDialogSettings();
    auto [new_atom, maybe_enhanced_stereo] = create_atom_with_properties(props);
    m_mol_model->mutateAtoms({m_atom}, *new_atom, maybe_enhanced_stereo);
    QDialog::accept();
}

/**
 * Set the atom property value to the value in the widget if and only if the
 * widget is enabled and not hidden.
 */
template <typename T> static void
set_to_dialog_value_if_widget_enabled(std::optional<T>& value,
                                      const BlankableSpinBox* const sb)
{
    // note that !widget->isHidden() is *not* the same as widget->isVisible().
    // The isVisible() method takes into account the visibility of
    // parent/grandparent/etc widgets, but isHidden() only considers the
    // visibility setting of the widget itself.
    if (sb->isEnabled() && !sb->isHidden()) {
        value = sb->optionalValue();
    }
}

template <typename T> static void
set_to_dialog_value_if_widget_enabled(T& value, const QSpinBox* const sb)
{
    if (sb->isEnabled() && !sb->isHidden()) {
        value = sb->value();
    }
}

template <typename T> static void
set_to_dialog_value_if_widget_enabled(T& value, const QComboBox* const combo)
{
    if (combo->isEnabled() && !combo->isHidden()) {
        value = combo->currentData().value<T>();
    }
}

std::shared_ptr<AbstractAtomProperties>
EditAtomPropertiesDialog::getDialogSettings() const
{
    if (ui->set_as_atom_rb->isChecked()) {
        AtomProperties atom_props;
        auto element_text = ui->atom_element_le->text();
        atom_props.element = get_element_from_text(element_text);

        set_to_dialog_value_if_widget_enabled(atom_props.isotope,
                                              ui->atom_isotope_sb);
        set_to_dialog_value_if_widget_enabled(atom_props.charge,
                                              ui->atom_charge_sb);
        set_to_dialog_value_if_widget_enabled(atom_props.unpaired_electrons,
                                              ui->atom_unpaired_sb);
        set_enhanced_stereo_properties(ui->atom_stereo_combo,
                                       ui->atom_stereo_sb, atom_props);
        return std::make_shared<AtomProperties>(atom_props);
    } else {

        AtomQueryProperties query_props;
        query_props.query_type =
            ui->query_type_combo->currentData().value<QueryType>();
        switch (query_props.query_type) {
            case QueryType::ALLOWED_LIST:
            case QueryType::NOT_ALLOWED_LIST: {
                auto element_list = ui->element_list_le->text();
                query_props.allowed_list = get_elements_from_list(element_list);
                break;
            }
            case QueryType::WILDCARD:
                query_props.wildcard =
                    ui->wildcard_combo->currentData().value<AtomQuery>();
                break;
            case QueryType::SPECIFIC_ELEMENT: {
                auto element_text = ui->query_element_le->text();
                query_props.element = get_element_from_text(element_text);
                break;
            }
            case QueryType::RGROUP:
                query_props.r_group = ui->rgroup_sb->value();
                break;
            case QueryType::SMARTS:
                query_props.smarts_query = getSmartsQuery();
                break;
        }
        set_to_dialog_value_if_widget_enabled(query_props.isotope,
                                              ui->query_isotope_sb);
        set_to_dialog_value_if_widget_enabled(query_props.charge,
                                              ui->query_charge_sb);
        set_to_dialog_value_if_widget_enabled(query_props.unpaired_electrons,
                                              ui->query_unpaired_sb);
        set_enhanced_stereo_properties(ui->query_stereo_combo,
                                       ui->query_stereo_sb, query_props);
        set_to_dialog_value_if_widget_enabled(query_props.total_h_type,
                                              ui->total_h_combo);
        set_to_dialog_value_if_widget_enabled(query_props.total_h_exact_val,
                                              ui->total_h_sb);
        set_to_dialog_value_if_widget_enabled(query_props.num_connections,
                                              ui->num_connections_sb);
        set_to_dialog_value_if_widget_enabled(query_props.aromaticity,
                                              ui->aromaticity_combo);
        set_to_dialog_value_if_widget_enabled(query_props.ring_count_type,
                                              ui->ring_count_combo);
        set_to_dialog_value_if_widget_enabled(query_props.ring_count_exact_val,
                                              ui->ring_count_sb);
        set_to_dialog_value_if_widget_enabled(query_props.ring_bond_count_type,
                                              ui->ring_bond_count_combo);
        set_to_dialog_value_if_widget_enabled(
            query_props.ring_bond_count_exact_val, ui->ring_bond_count_sb);
        set_to_dialog_value_if_widget_enabled(query_props.smallest_ring_size,
                                              ui->smallest_ring_size_sb);
        return std::make_shared<AtomQueryProperties>(query_props);
    }
}

std::string EditAtomPropertiesDialog::getSmartsQuery() const
{
    auto smarts = ui->smarts_query_le->text().trimmed().toStdString();
    if (smarts.empty()) {
        throw InvalidAtomPropertyError("Missing SMARTS query.");
    }
    auto* atom = RDKit::SmartsToAtom(smarts);
    if (atom == nullptr) {
        throw InvalidAtomPropertyError("Invalid SMARTS query.");
    }
    delete atom;
    return smarts;
}

void EditAtomPropertiesDialog::loadProperties(
    const std::shared_ptr<AbstractAtomProperties> props)
{
    if (props->isQuery()) {
        ui->set_as_query_rb->click();
        loadAtomProperties(AtomProperties());
        // reset the query page to defaults
        loadQueryProperties(
            *static_cast<const AtomQueryProperties*>(props.get()));
    } else {
        ui->set_as_atom_rb->click();
        // reset the atom page to defaults
        loadAtomProperties(*static_cast<const AtomProperties*>(props.get()));
        loadQueryProperties(AtomQueryProperties());
    }

    // load the enhanced stereo properties into both pages so that the label and
    // combo box are correctly enabled or disabled
    set_enhanced_stereo_widgets(ui->atom_stereo_lbl, ui->atom_stereo_combo,
                                ui->atom_stereo_sb, *props);
    set_enhanced_stereo_widgets(ui->query_stereo_lbl, ui->query_stereo_combo,
                                ui->query_stereo_sb, *props);

    updateButtonsEnabled();
}

void EditAtomPropertiesDialog::loadAtomProperties(
    const AtomProperties& atom_props)
{
    auto atomic_symbol = get_atomic_symbol_from_element(atom_props.element);
    ui->atom_element_le->setText(atomic_symbol);
    ui->atom_isotope_sb->setOptionalValue(atom_props.isotope);
    ui->atom_charge_sb->setValue(atom_props.charge);
    ui->atom_unpaired_sb->setValue(atom_props.unpaired_electrons);
}

void EditAtomPropertiesDialog::loadQueryProperties(
    const AtomQueryProperties& query_props)
{
    set_combo_box_data(ui->query_type_combo, query_props.query_type);
    clearQueryTypeFields();
    switch (query_props.query_type) {
        case QueryType::ALLOWED_LIST:
        case QueryType::NOT_ALLOWED_LIST: {
            auto allowed_list =
                join_all_atomic_symbols(query_props.allowed_list, ", ");
            ui->element_list_le->setText(allowed_list);
            break;
        }
        case QueryType::WILDCARD:
            set_combo_box_data(ui->wildcard_combo, query_props.wildcard);
            break;
        case QueryType::SPECIFIC_ELEMENT: {
            auto atomic_symbol =
                get_atomic_symbol_from_element(query_props.element);
            ui->query_element_le->setText(atomic_symbol);
            break;
        }
        case QueryType::RGROUP:
            ui->rgroup_sb->setValue(query_props.r_group);
            break;
        case QueryType::SMARTS: {
            auto qsmarts = QString::fromStdString(query_props.smarts_query);
            ui->smarts_query_le->setText(qsmarts);
            break;
        }
    }

    ui->query_isotope_sb->setOptionalValue(query_props.isotope);
    ui->query_charge_sb->setOptionalValue(query_props.charge);
    ui->query_unpaired_sb->setOptionalValue(query_props.unpaired_electrons);

    set_combo_box_data(ui->total_h_combo, query_props.total_h_type);
    ui->total_h_sb->setValue(query_props.total_h_exact_val);
    ui->num_connections_sb->setOptionalValue(query_props.num_connections);
    set_combo_box_data(ui->aromaticity_combo, query_props.aromaticity);
    set_combo_box_data(ui->ring_count_combo, query_props.ring_count_type);
    ui->ring_count_sb->setValue(query_props.ring_count_exact_val);
    set_combo_box_data(ui->ring_bond_count_combo,
                       query_props.ring_bond_count_type);
    ui->ring_bond_count_sb->setValue(query_props.ring_bond_count_exact_val);
    ui->smallest_ring_size_sb->setOptionalValue(query_props.smallest_ring_size);
}

void EditAtomPropertiesDialog::clearQueryTypeFields()
{
    ui->element_list_le->clear();
    ui->query_element_le->clear();
    ui->rgroup_sb->setValue(1);
    ui->wildcard_combo->setCurrentIndex(0);
}

/**
 * Transfer enhanced stereo settings between the "Set as: Atom" page and the
 * "Set as: Query" page.
 */
static void transfer_enhanced_stereo_settings(const QComboBox* const from_combo,
                                              const QSpinBox* const from_sb,
                                              QComboBox* const to_combo,
                                              QSpinBox* const to_sb)
{
    if (to_combo->isEnabled()) {
        to_combo->setCurrentIndex(from_combo->currentIndex());
        to_sb->setValue(from_sb->value());
    }
}

/**
 * Transfer settings from a QSpinBox to a BlankableSpinBox. This method will
 * not overwrite a blank value with a zero since there's no way to distinguish
 * those two settings in a non-blankable spin box. In that scenario, this
 * method will be a no-op.
 */
static void transfer_sb_to_optional_sb(const QSpinBox* const from_sb,
                                       BlankableSpinBox* const to_sb)
{
    auto value = from_sb->value();
    if (value != 0 || to_sb->optionalValue().has_value()) {
        to_sb->setValue(value);
    }
}

void EditAtomPropertiesDialog::onAtomOrQueryToggled(const int page_id)
{
    auto page = ui->atom_query_stacked_wdg->widget(page_id);
    auto query_type = ui->query_type_combo->currentData().value<QueryType>();
    if (page == ui->edit_atom_page) {
        // we're switching to the atom page
        if (query_type == QueryType::SPECIFIC_ELEMENT) {
            ui->atom_element_le->setText(ui->query_element_le->text());
        } else if (query_type == QueryType::ALLOWED_LIST ||
                   query_type == QueryType::NOT_ALLOWED_LIST) {
            transferElementListToSpecificElement(ui->atom_element_le);
        }
        ui->atom_isotope_sb->setOptionalValue(
            ui->query_isotope_sb->optionalValue());
        ui->atom_charge_sb->setValue(
            ui->query_charge_sb->optionalValue().value_or(0));
        ui->atom_unpaired_sb->setValue(
            ui->query_unpaired_sb->optionalValue().value_or(0));
        transfer_enhanced_stereo_settings(
            ui->query_stereo_combo, ui->query_stereo_sb, ui->atom_stereo_combo,
            ui->atom_stereo_sb);
    } else {
        // we're switching to the query page
        if (query_type == QueryType::SPECIFIC_ELEMENT) {
            auto atom_element_text = ui->atom_element_le->text();
            try {
                get_element_from_text(atom_element_text);
                ui->query_element_le->setText(atom_element_text);
            } catch (InvalidAtomPropertyError&) {
                // the element isn't valid, so don't bother to transfer it
            }
        } else if (query_type == QueryType::ALLOWED_LIST ||
                   query_type == QueryType::NOT_ALLOWED_LIST) {
            transferSpecificElementToElementList(ui->atom_element_le);
        }
        ui->query_isotope_sb->setOptionalValue(
            ui->atom_isotope_sb->optionalValue());
        transfer_sb_to_optional_sb(ui->atom_charge_sb, ui->query_charge_sb);
        transfer_sb_to_optional_sb(ui->atom_unpaired_sb, ui->query_unpaired_sb);
        transfer_enhanced_stereo_settings(
            ui->atom_stereo_combo, ui->atom_stereo_sb, ui->query_stereo_combo,
            ui->query_stereo_sb);
    }
    updateButtonsEnabled();
}

void EditAtomPropertiesDialog::transferSpecificElementToElementList(
    const QLineEdit* const element_le)
{
    auto element_text = element_le->text();
    try {
        get_element_from_text(element_text);
    } catch (InvalidAtomPropertyError&) {
        // the line edit doesn't specify a valid element, so don't bother to
        // transfer it
        return;
    }
    // don't transfer the element list if the element list already starts with
    // that element, in case the specific element field was populated by
    // truncating the element list
    auto element_list = split_element_list(ui->element_list_le->text());
    if (element_list.isEmpty() || element_text != element_list[0]) {
        ui->element_list_le->setText(element_text);
    }
}

void EditAtomPropertiesDialog::transferElementListToSpecificElement(
    QLineEdit* const element_le)
{
    auto element_list = split_element_list(ui->element_list_le->text());
    if (element_list.isEmpty()) {
        // there's nothing to transfer
        return;
    }
    auto first_atomic_symbol = element_list[0];
    try {
        get_element_from_text(first_atomic_symbol);
    } catch (InvalidAtomPropertyError&) {
        // the atomic symbol is invalid, so don't bother to transfer it
        return;
    }
    element_le->setText(first_atomic_symbol);
}

void EditAtomPropertiesDialog::updateButtonsEnabled()
{
    // if the timer was started and then this method was called explicitly,
    // there's no need to call it again when the timer goes off
    m_update_buttons_enabled_timer->stop();

    bool enable_ok = true;
    bool enable_reset = true;
    std::string ok_tooltip;
    std::string reset_tooltip;
    std::shared_ptr<AbstractAtomProperties> props = nullptr;

    // check if the settings are valid
    try {
        props = getDialogSettings();
    } catch (InvalidAtomPropertyError& e) {
        enable_ok = false;
        ok_tooltip = e.what();
    }

    // check if the user has changed any settings
    if (enable_ok && *props == *read_properties_from_atom(m_atom)) {
        enable_ok = false;
        enable_reset = false;
        ok_tooltip = reset_tooltip = "No settings have been changed.";
    }

    auto* ok_btn = ui->buttonBox->button(QDialogButtonBox::Ok);
    ok_btn->setEnabled(enable_ok);
    ok_btn->setToolTip(QString::fromStdString(ok_tooltip));

    auto* reset_btn = ui->buttonBox->button(QDialogButtonBox::Reset);
    reset_btn->setEnabled(enable_reset);
    reset_btn->setToolTip(QString::fromStdString(reset_tooltip));
}

void EditAtomPropertiesDialog::onQueryTypeComboBoxChanged(const int combo_index)
{
    auto variant = ui->query_type_combo->itemData(combo_index);
    QWidget* switch_to_page = nullptr;
    bool show_periodic_table_btn = true;
    bool enable_isotope = true;
    bool enable_charge = true;
    bool enable_unpaired_electrons = true;
    bool enable_advanced_tab = true;
    switch (variant.value<QueryType>()) {
        case QueryType::ALLOWED_LIST:
        case QueryType::NOT_ALLOWED_LIST:
            switch_to_page = ui->element_list_page;
            enable_isotope = false;
            enable_unpaired_electrons = false;
            break;
        case QueryType::WILDCARD:
            switch_to_page = ui->wildcard_page;
            show_periodic_table_btn = false;
            enable_unpaired_electrons = false;
            break;
        case QueryType::SPECIFIC_ELEMENT:
            switch_to_page = ui->query_element_page;
            break;
        case QueryType::RGROUP:
            switch_to_page = ui->rgroup_page;
            show_periodic_table_btn = false;
            enable_isotope = false;
            enable_charge = false;
            enable_unpaired_electrons = false;
            break;
        case QueryType::SMARTS:
            switch_to_page = ui->smarts_query_page;
            show_periodic_table_btn = false;
            enable_isotope = false;
            enable_charge = false;
            enable_unpaired_electrons = false;
            enable_advanced_tab = false;
            break;
    }

    // transfer elements between specific element and element list
    auto switch_from_page = ui->query_type_stacked_widget->currentWidget();
    if (switch_from_page == ui->query_element_page) {
        transferSpecificElementToElementList(ui->query_element_le);
    } else if (switch_from_page == ui->element_list_page &&
               switch_to_page != ui->element_list_page) {
        transferElementListToSpecificElement(ui->query_element_le);
    }

    ui->query_type_stacked_widget->setCurrentWidget(switch_to_page);
    ui->query_periodic_table_btn->setVisible(show_periodic_table_btn);
    ui->query_isotope_lbl->setEnabled(enable_isotope);
    ui->query_isotope_sb->setEnabled(enable_isotope);
    ui->query_charge_lbl->setEnabled(enable_charge);
    ui->query_charge_sb->setEnabled(enable_charge);
    ui->query_unpaired_lbl->setEnabled(enable_unpaired_electrons);
    ui->query_unpaired_sb->setEnabled(enable_unpaired_electrons);
    ui->advanced_tab->setEnabled(enable_advanced_tab);
    auto advanced_tab_idx = ui->edit_query_tab_wdg->indexOf(ui->advanced_tab);
    ui->edit_query_tab_wdg->setTabEnabled(advanced_tab_idx,
                                          enable_advanced_tab);
    updateButtonsEnabled();
}

void EditAtomPropertiesDialog::onAtomPeriodicTableElementSelected(
    Element element)
{
    auto symbol = get_atomic_symbol_from_element(element);
    ui->atom_element_le->setText(symbol);
    ui->atom_periodic_table_btn->getPopupWidget()->close();
}

void EditAtomPropertiesDialog::onQueryPeriodicTableElementSelected(
    Element element)
{
    auto symbol = get_atomic_symbol_from_element(element);
    auto query_type = ui->query_type_combo->currentData().value<QueryType>();
    if (query_type == QueryType::SPECIFIC_ELEMENT) {
        ui->query_element_le->setText(symbol);
        ui->query_periodic_table_btn->getPopupWidget()->close();
    } else {
        auto text = ui->element_list_le->text();
        // trim unnecessary whitespace
        text = text.simplified();
        if (!text.isEmpty() && !text.endsWith(',')) {
            text.append(", ");
        }
        text.append(symbol);
        ui->element_list_le->setText(text);
    }
#ifdef __EMSCRIPTEN__
    // On WASM builds, the dialog will incorrectly be moved in front of the
    // periodic table whenever the dialog's OK button gets re-enabled, so make
    // sure that we put the periodic table back in front
    ui->query_periodic_table_btn->getPopupWidget()->raise();
#endif
}

ElementValidator::ElementValidator(const bool allow_list, QObject* parent) :
    QValidator(parent),
    m_allow_list(allow_list)
{
}

QValidator::State ElementValidator::validate(QString& input, int& pos) const
{
    // check for invalid characters
    for (auto& cur_char : input) {
        if (!cur_char.isLetter() &&
            !(m_allow_list && (cur_char == ',' || cur_char == ' '))) {
            // don't allow the invalid characters to be added to the line edit
            return QValidator::Invalid;
        }
    }

    // find out if the atomic symbols are valid (i.e. are they real elements)
    try {
        get_elements_from_list(input);
    } catch (InvalidAtomPropertyError&) {
        return QValidator::Intermediate;
    }
    return QValidator::Acceptable;
}

} // namespace sketcher
} // namespace schrodinger
