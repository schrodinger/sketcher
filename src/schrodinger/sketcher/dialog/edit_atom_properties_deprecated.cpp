#include "schrodinger/sketcher/dialog/edit_atom_properties_deprecated.h"

#include <set>
#include <sstream>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <QPushButton>

#include "schrodinger/sketcher/Atom.h"
#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/atom_type_mappings.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/ui/ui_common_atom_properties_widget.h"
#include "schrodinger/sketcher/ui/ui_edit_atom_properties_deprecated.h"
#include "schrodinger/sketcher/undoable_structure_change_wrapper.h"
#include "schrodinger/sketcher/widget/periodic_table_widget.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::QueryType);
Q_DECLARE_METATYPE(schrodinger::sketcher::AtomQuery);

namespace
{

using schrodinger::sketcher::AtomQuery;
using schrodinger::sketcher::QueryType;

const QString ANY_COMBO_KEY{"(any)"};
const QString AROMATIC_COMBO_KEY{"Aromatic (a)"};
const QString ALIPHATIC_COMBO_KEY{"Aliphatic (A)"};
const QString NONZERO_COMBO_KEY{">0"};
const QString EXACT_COMBO_KEY{"exactly"};

const std::vector<std::pair<QueryType, QString>> query_type_titles = {
    {QueryType::ALLOWED_LIST, "Allowed List"},
    {QueryType::NOT_ALLOWED_LIST, "Not Allowed List"},
    {QueryType::WILDCARD, "Wildcard"},
    {QueryType::SPECIFIC_ELEMENT, "Specific Element"},
    {QueryType::RGROUP, "R-Group"},
};

const std::vector<std::pair<AtomQuery, QString>> atomquery_titles = {
    {AtomQuery::A, "Any heavy atom (A)"}, {AtomQuery::Q, "Heteroatom (Q)"},
    {AtomQuery::M, "Metal (M)"},          {AtomQuery::X, "Halogen (X)"},
    {AtomQuery::AH, "Any or H (AH)"},     {AtomQuery::QH, "Hetero or H (QH)"},
    {AtomQuery::MH, "Metal or H (MH)"},   {AtomQuery::XH, "Halogen or H (XH)"},
};

void update_ring_count_spinbox(QSpinBox* sb, const QString& combo_key_text)
{
    if (combo_key_text == EXACT_COMBO_KEY) {
        sb->setValue(sb->value());
        sb->setEnabled(true);
    } else {
        sb->clear();
        sb->setEnabled(false);
    }
}

bool has_atom_specified_chirality(const RDKit::Atom* atom)
{
    auto sinfo = RDKit::Chirality::detail::getStereoInfo(atom);
    return sinfo.type == RDKit::Chirality::StereoType::Atom_Tetrahedral &&
           sinfo.specified == RDKit::Chirality::StereoSpecified::Specified;
}

} // anonymous namespace

namespace schrodinger
{
namespace sketcher
{

const uint SET_AS_ATOM_ID = 0;
const uint SET_AS_QUERY_ID = 1;

CommonAtomPropertiesWidget::CommonAtomPropertiesWidget(QWidget* parent) :
    QWidget(parent)
{
    ui.reset(new Ui::CommonAtomPropertiesWidget());
    ui->setupUi(this);

    const std::map<QString, ::EnhancedStereoType> items = {
        {"ABS", ::EnhancedStereoType::ABS},
        {"AND", ::EnhancedStereoType::AND},
        {"OR", ::EnhancedStereoType::OR},
    };
    for (const auto& item_pair : items) {
        ui->enhanced_stereo_combo->addItem(
            item_pair.first, QVariant::fromValue(item_pair.second));
    }

    connect(ui->enhanced_stereo_combo,
            QOverload<int>::of(&QComboBox::currentIndexChanged), this,
            &CommonAtomPropertiesWidget::updateEnhancedStereoSpinBox);

    setEnhancedStereoTypeComboValue(::EnhancedStereoType::ABS);

    // Cap charges
    ui->charge_sb->setMinimum(-ATOM_CHARGE_LIMIT);
    ui->charge_sb->setMaximum(ATOM_CHARGE_LIMIT);
}

void CommonAtomPropertiesWidget::setEnhancedStereoTypeComboValue(
    const ::EnhancedStereoType& enhanced_stereo_type)
{
    auto index = ui->enhanced_stereo_combo->findData(
        QVariant::fromValue(enhanced_stereo_type));
    ui->enhanced_stereo_combo->setCurrentIndex(index);
    updateEnhancedStereoSpinBox();
}

void CommonAtomPropertiesWidget::onQueryTypeComboValueChanged(
    QueryType query_type)
{
    auto in_specific_element = query_type == QueryType::SPECIFIC_ELEMENT;
    auto in_wildcard = query_type == QueryType::WILDCARD;

    bool isotope_enabled = in_specific_element || in_wildcard;
    ui->isotope_lbl->setEnabled(isotope_enabled);
    ui->isotope_le->setEnabled(isotope_enabled);

    bool charge_enabled = in_specific_element || in_wildcard ||
                          query_type == QueryType::ALLOWED_LIST ||
                          query_type == QueryType::NOT_ALLOWED_LIST;
    ui->charge_lbl->setEnabled(charge_enabled);
    ui->charge_sb->setEnabled(charge_enabled);

    ui->unpaired_electrons_lbl->setEnabled(in_specific_element);
    ui->unpaired_elec_sb->setEnabled(in_specific_element);
}

void CommonAtomPropertiesWidget::setEnhancedStereoTypeComboEnabled(bool enable)
{
    ui->enhanced_stereo_lbl->setEnabled(enable);
    ui->enhanced_stereo_combo->setEnabled(enable);
}

CommonAtomPropertiesWidget::~CommonAtomPropertiesWidget() = default;

void CommonAtomPropertiesWidget::readAtomInfo(const sketcherAtom& atom)
{
    ui->isotope_le->setText(QString::number(atom.getIsotope()));
    if (atom.getIsotope() == 0) {
        ui->isotope_le->clear();
    }
    ui->charge_sb->setValue(atom.getCharge());
    ui->unpaired_elec_sb->setValue(atom.getUnpairedElectronsN());

    // Atoms must have a defined chirality in order to get enhanced stereo
    auto rdk_atom = atom.getRDKAtom();
    auto enhanced_stereo_type = ::EnhancedStereoType::ABS;
    auto enhanced_stereo_group_id = 0;
    if (rdk_atom != nullptr && has_atom_specified_chirality(rdk_atom)) {
        auto enh_stereo = atom.getEnhancedStereo();
        enhanced_stereo_type = enh_stereo.type;
        enhanced_stereo_group_id = enh_stereo.group_id;
    }

    // Set the combo value last because altering this combo box will signal
    // `updateEnhancedStereoSpinBox()` which may then update the enhanced stereo
    // spin box anyway.
    ui->enhanced_stereo_sb->setValue(enhanced_stereo_group_id);
    setEnhancedStereoTypeComboValue(enhanced_stereo_type);
}

void CommonAtomPropertiesWidget::writeAtomInfo(sketcherAtom& atom) const
{
    auto get_value = [](auto widget) {
        if (!widget->isEnabled() || widget->text().isEmpty()) {
            return 0; // Reset if the widget is inappropriate (disabled)
        }
        return widget->text().toInt();
    };

    // TODO: Add isotope validator though RDKit::getAbundanceForIsotope
    int ui_isotope_value = get_value(ui->isotope_le);

    atom.setIsotope(std::min(std::abs(ui_isotope_value), MAX_ISOTOPE_VALUE));
    atom.setCharge(get_value(ui->charge_sb));
    atom.setUnpairedElectronsN(get_value(ui->unpaired_elec_sb));

    auto type =
        ui->enhanced_stereo_combo->currentData().value<::EnhancedStereoType>();
    unsigned int group_id = ui->enhanced_stereo_sb->value();
    if (type == ::EnhancedStereoType::ABS) {
        group_id = 0;
    }
    atom.setEnhancedStereo({type, group_id});
}

void CommonAtomPropertiesWidget::updateEnhancedStereoSpinBox()
{
    auto enhanced_stereo_type =
        ui->enhanced_stereo_combo->currentData().value<::EnhancedStereoType>();
    switch (enhanced_stereo_type) {
        // Absolute enhanced stereo does not allow groups
        case ::EnhancedStereoType::ABS:
            ui->enhanced_stereo_sb->clear();
            ui->enhanced_stereo_sb->setEnabled(false);
            break;
        case ::EnhancedStereoType::AND:
        case ::EnhancedStereoType::OR:
            // show value should it have been previously cleared
            ui->enhanced_stereo_sb->setValue(ui->enhanced_stereo_sb->value());
            ui->enhanced_stereo_sb->setEnabled(true);
            break;
    }
}

template <typename Func1, typename Func2> void
link_widgets(const typename QtPrivate::FunctionPointer<Func1>::Object* widget_a,
             const typename QtPrivate::FunctionPointer<Func2>::Object* widget_b,
             Func1 signal, Func2 slot)
{
    QObject::connect(widget_a, signal, widget_b, slot);
    QObject::connect(widget_b, signal, widget_a, slot);
}

EditAtomPropertiesDialogDeprecated::EditAtomPropertiesDialogDeprecated(
    SketcherModel* model, sketcherAtom& atom, QWidget* parent,
    Qt::WindowFlags f) :
    ModalDialog(parent, f),
    m_atom(atom),
    m_sketcher_model(model)
{
    if (model == nullptr) {
        throw std::runtime_error("This dialog cannot be created without a"
                                 " model");
    }

    ui.reset(new Ui::EditAtomPropertiesDialogDeprecated());
    setupDialogUI(*ui);

    ui->set_as_group->setId(ui->set_as_atom_rb, SET_AS_ATOM_ID);
    ui->set_as_group->setId(ui->set_as_query_rb, SET_AS_QUERY_ID);

    // Connect radio buttons to stacked widget
    auto stack_wgt = ui->atom_query_stacked_wdg;
    connect(ui->set_as_atom_rb, &QRadioButton::clicked, stack_wgt,
            std::bind(&QStackedWidget::setCurrentWidget, stack_wgt,
                      ui->edit_atom_page));
    connect(ui->set_as_query_rb, &QRadioButton::clicked, stack_wgt,
            std::bind(&QStackedWidget::setCurrentWidget, stack_wgt,
                      ui->edit_query_page));

    auto resetButton = ui->buttonBox->button(QDialogButtonBox::Reset);
    connect(resetButton, &QPushButton::clicked, this,
            &EditAtomPropertiesDialogDeprecated::reset);

    auto& atom_ui = ui->atom_common_props_wdg->ui;
    auto& query_ui = ui->query_common_props_wdg->ui;

    // Link the properties common to both the atom and query pages
    link_widgets(atom_ui->isotope_le, query_ui->isotope_le,
                 &QLineEdit::textChanged, &QLineEdit::setText);

    link_widgets(atom_ui->enhanced_stereo_combo,
                 query_ui->enhanced_stereo_combo,
                 QOverload<int>::of(&QComboBox::currentIndexChanged),
                 &QComboBox::setCurrentIndex);

    const std::vector<std::pair<QSpinBox*, QSpinBox*>> spin_boxes = {
        {atom_ui->charge_sb, query_ui->charge_sb},
        {atom_ui->unpaired_elec_sb, query_ui->unpaired_elec_sb},
        {atom_ui->enhanced_stereo_sb, query_ui->enhanced_stereo_sb}};

    for (const auto& sb : spin_boxes) {
        link_widgets(sb.first, sb.second,
                     QOverload<int>::of(&QSpinBox::valueChanged),
                     &QSpinBox::setValue);
    }

    // Set up the periodic table button/widget
    m_periodic_table_wdg = new PeriodicTableWidget(this);
    m_periodic_table_wdg->setWindowFlags(Qt::Popup);
    connect(
        m_periodic_table_wdg, &PeriodicTableWidget::elementSelected, this,
        &EditAtomPropertiesDialogDeprecated::onPeriodicTableElementSelected);
    ui->periodic_table_btn->setPopupDelay(0);
    ui->periodic_table_btn->setPopupWidget(m_periodic_table_wdg);
    ui->periodic_table_btn->showPopupIndicator(false);

    // Set up the query type combo box
    for (auto& pair : ::query_type_titles) {
        ui->query_type_combo->addItem(pair.second,
                                      QVariant::fromValue(pair.first));
    }
    connect(ui->query_type_combo,
            QOverload<int>::of(&QComboBox::currentIndexChanged), this,
            &EditAtomPropertiesDialogDeprecated::onQueryTypeComboValueChanged);
    setQueryTypeComboValue(QueryType::SPECIFIC_ELEMENT); // Default
    onQueryTypeComboValueChanged();

    // Set up the wildcard combo box
    for (auto& pair : ::atomquery_titles) {
        ui->wildcard_combo->addItem(pair.second,
                                    QVariant::fromValue(pair.first));
    }

    auto rdk_atom = m_atom.getRDKAtom();
    bool enable = rdk_atom != nullptr && has_atom_specified_chirality(rdk_atom);
    ui->atom_common_props_wdg->setEnhancedStereoTypeComboEnabled(enable);
    ui->query_common_props_wdg->setEnhancedStereoTypeComboEnabled(enable);

    connect(ui->ring_count_combo, &QComboBox::currentTextChanged, this,
            std::bind(::update_ring_count_spinbox, ui->ring_count_sb,
                      std::placeholders::_1));
    connect(ui->ring_bond_count_combo, &QComboBox::currentTextChanged, this,
            std::bind(::update_ring_count_spinbox, ui->ring_bond_count_sb,
                      std::placeholders::_1));

    // Update whether the OK button is enabled when relevant widgets change
    std::vector<QLineEdit*> objects{ui->specific_element_le,
                                    ui->element_list_le, ui->element_le};
    for (auto object : objects) {
        connect(object, &QLineEdit::textChanged, this,
                &EditAtomPropertiesDialogDeprecated::updateOKButtonEnabled);
    }

    connect(ui->set_as_atom_rb, &QRadioButton::toggled, this,
            &EditAtomPropertiesDialogDeprecated::updateOKButtonEnabled);
    connect(ui->set_as_query_rb, &QRadioButton::toggled, this,
            &EditAtomPropertiesDialogDeprecated::updateOKButtonEnabled);

    reset();
}

EditAtomPropertiesDialogDeprecated::~EditAtomPropertiesDialogDeprecated() =
    default;

void EditAtomPropertiesDialogDeprecated::accept()
{
    UndoableStructureChangeWrapper undo_wrapper(*m_atom.getSketcherScene(),
                                                "Edit Atom Properties");
    auto current_page = ui->atom_query_stacked_wdg->currentWidget();

    try {
        if (current_page == ui->edit_atom_page) {
            writeAtomInfo();
        } else { // ui->edit_query_page
            writeQueryInfo();
        }
    } catch (const std::runtime_error&) {
        // The provided element text is not understood
    }
    m_atom.getSketcherScene()->invalidateStructure();
    QDialog::accept();
}

void EditAtomPropertiesDialogDeprecated::reset()
{
    // Reset the "Set as" radio buttons
    if (m_atom.isQuery() || m_atom.hasAdvancedQueryFeatures()) {
        ui->set_as_query_rb->click();
    } else {
        ui->set_as_atom_rb->click();
    }

    // Reset the "Set as atom" subwidgets
    auto atom_type = m_atom.getAtomType();
    std::string element_text;
    if (is_atomic_number(atom_type)) {
        element_text = atomic_number_to_symbol(atom_type);
    }
    ui->element_le->setText(QString::fromStdString(element_text));

    ui->atom_common_props_wdg->readAtomInfo(m_atom);

    resetGeneralQuerySubwidgets();
    resetAdvancedQuerySubwidgets();
}

void EditAtomPropertiesDialogDeprecated::resetGeneralQuerySubwidgets()
{
    auto query_type = QueryType::SPECIFIC_ELEMENT;
    if (m_atom.isAllowedAtomListQuery()) {
        query_type = QueryType::ALLOWED_LIST;
    } else if (m_atom.isDisallowedAtomListQuery()) {
        query_type = QueryType::NOT_ALLOWED_LIST;
    } else if (m_atom.isRGroup()) {
        query_type = QueryType::RGROUP;
    } else if (m_atom.isAnyAtomWildcard() || m_atom.isHeavyAtomWildcard() ||
               m_atom.isAnyNonCarbonWildcard() ||
               m_atom.isNonCarbonHeavyAtomWildcard() ||
               m_atom.isMetalOrHWildcard() || m_atom.isMetalWildcard() ||
               m_atom.isHalogenOrHWildcard() || m_atom.isHalogenWildcard()) {
        query_type = QueryType::WILDCARD;
    }
    setQueryTypeComboValue(query_type);

    auto atom_type = m_atom.getAtomType();
    auto wildcard_combo_idx = 0;
    if (query_type == QueryType::WILDCARD) {
        auto atom_type_qvariant = QVariant::fromValue(
            get_atomquery_for_atomtype(AtomTypes(atom_type)));
        wildcard_combo_idx = ui->wildcard_combo->findData(atom_type_qvariant);
    }
    ui->wildcard_combo->setCurrentIndex(wildcard_combo_idx);

    std::string element_list_text;
    std::string element_text;
    if (is_atomic_number(atom_type)) {
        element_text = atomic_number_to_symbol(atom_type);
    }
    if (!m_atom.isQuery()) {
        element_list_text = element_text;
    } else if (m_atom.isAllowedAtomListQuery() ||
               m_atom.isDisallowedAtomListQuery()) {
        element_list_text = m_atom.getAtomicNumbersListString();
    }
    ui->element_list_le->setText(QString::fromStdString(element_list_text));

    ui->specific_element_le->setText(QString::fromStdString(element_text));

    auto r_group_number = m_sketcher_model->getNextRGroupNumber();
    if (m_atom.isRGroup()) {
        r_group_number = m_atom.getRGroupNumber();
    }
    ui->rgroup_sb->setValue(r_group_number);

    ui->query_common_props_wdg->readAtomInfo(m_atom);
}

void EditAtomPropertiesDialogDeprecated::resetAdvancedQuerySubwidgets()
{
    std::string atom_smarts;
    m_atom.getStringProperty(ATOM_SMARTS_KEY, atom_smarts);

    QString total_h_text = "";
    if (atom_smarts.empty()) {
        int total_h = -1;
        if (m_atom.getIntProperty(TOTAL_H_COUNT_KEY, total_h) && total_h > -1) {
            total_h_text = QString::number(total_h);
        }
    }
    ui->total_H_le->setText(total_h_text);

    QString num_connections_text{""};
    if (atom_smarts.empty()) {
        int num_connections = -1;
        if (m_atom.getIntProperty(NUMBER_OF_CONNECTIONS_KEY, num_connections) &&
            num_connections > -1) {
            num_connections_text = QString::number(num_connections);
        }
    }
    ui->num_connections_le->setText(num_connections_text);

    auto aromaticity_key = ANY_COMBO_KEY;
    if (atom_smarts.empty()) {
        int aromaticity = -1;
        if (m_atom.getIntProperty(AROMATICITY_KEY, aromaticity)) {
            if (aromaticity == sketcherAtom::ALIPHATIC_ATOM) {
                aromaticity_key = ALIPHATIC_COMBO_KEY;
            } else if (aromaticity == sketcherAtom::AROMATIC_ATOM) {
                aromaticity_key = AROMATIC_COMBO_KEY;
            }
        }
    }
    auto aromaticity_index = ui->aromaticity_combo->findText(aromaticity_key);
    ui->aromaticity_combo->setCurrentIndex(aromaticity_index);

    auto ring_count_key = ANY_COMBO_KEY;
    int ring_count = 0;
    if (atom_smarts.empty()) {
        int min_rings = -1;
        if (m_atom.getIntProperty(MINIMUM_NUMBER_OF_RINGS_KEY, min_rings)) {
            int max_rings = -1;
            if (m_atom.getIntProperty(MAXIMUM_NUMBER_OF_RINGS_KEY, max_rings)) {
                if (min_rings == max_rings && min_rings > -1) {
                    ring_count_key = EXACT_COMBO_KEY;
                    ring_count = min_rings;
                }
            } else if (min_rings == 1 && max_rings == -1) {
                ring_count_key = NONZERO_COMBO_KEY;
            }
        }
    }
    ui->ring_count_sb->setValue(ring_count);
    auto ring_count_index = ui->ring_count_combo->findText(ring_count_key);
    ui->ring_count_combo->setCurrentIndex(ring_count_index);
    update_ring_count_spinbox(ui->ring_count_sb, ring_count_key);

    auto ring_bonds_key = ANY_COMBO_KEY;
    int ring_bond_count = 0;
    if (atom_smarts.empty()) {
        int ring_bonds = -1;
        if (m_atom.getIntProperty(RING_BOND_COUNT_KEY, ring_bonds) &&
            ring_bonds > -1) {
            ring_bonds_key = EXACT_COMBO_KEY;
            ring_bond_count = ring_bonds;
        }
    }
    ui->ring_bond_count_sb->setValue(ring_bond_count);
    auto ring_bonds_index = ui->ring_bond_count_combo->findText(ring_bonds_key);
    ui->ring_bond_count_combo->setCurrentIndex(ring_bonds_index);
    update_ring_count_spinbox(ui->ring_bond_count_sb, ring_bonds_key);

    QString ring_size_text{""};
    if (atom_smarts.empty()) {
        int ring_size = -1;
        if (m_atom.getIntProperty(MINIMUM_SIZE_OF_RINGS_KEY, ring_size) &&
            ring_size > -1) {
            ring_size_text = QString::number(ring_size);
        }
    }
    ui->smallest_ring_size_le->setText(ring_size_text);

    ui->smarts_query_le->setText(QString::fromStdString(atom_smarts));
}

QueryType EditAtomPropertiesDialogDeprecated::getQueryTypeComboValue() const
{
    return ui->query_type_combo->currentData().value<QueryType>();
}

AtomQuery EditAtomPropertiesDialogDeprecated::getWildcardComboValue() const
{
    return ui->wildcard_combo->currentData().value<AtomQuery>();
}

void EditAtomPropertiesDialogDeprecated::setQueryTypeComboValue(QueryType etype)
{
    auto idx = ui->query_type_combo->findData(QVariant::fromValue(etype));
    ui->query_type_combo->setCurrentIndex(idx);
}

void EditAtomPropertiesDialogDeprecated::onQueryTypeComboValueChanged()
{
    auto query_type = getQueryTypeComboValue();
    auto in_specific_element = query_type == QueryType::SPECIFIC_ELEMENT;
    auto in_element_list = query_type == QueryType::ALLOWED_LIST ||
                           query_type == QueryType::NOT_ALLOWED_LIST;

    ui->element_list_le->setVisible(in_element_list);
    ui->specific_element_le->setVisible(in_specific_element);
    ui->wildcard_combo->setVisible(query_type == QueryType::WILDCARD);
    ui->rgroup_sb->setVisible(query_type == QueryType::RGROUP);

    // Only show the periodic table when applicable
    m_periodic_table_wdg->setCloseOnClick(in_specific_element);
    ui->periodic_table_btn->setVisible(in_specific_element || in_element_list);

    // Only certain items are sensible for each query type
    ui->query_common_props_wdg->onQueryTypeComboValueChanged(query_type);
    updateOKButtonEnabled();
}

void EditAtomPropertiesDialogDeprecated::onPeriodicTableElementSelected(
    Element element)
{
    auto symbol = QString::fromStdString(
        atomic_number_to_symbol(static_cast<int>(element)));
    if (getQueryTypeComboValue() == QueryType::SPECIFIC_ELEMENT) {
        // Replace text with element
        ui->specific_element_le->setText(symbol);
    } else {
        QString text = ui->element_list_le->text();
        // Append element symbol to existing text
        if (!text.isEmpty()) {
            text += ", ";
        }
        text += symbol;
        ui->element_list_le->setText(text);
    }
}

void EditAtomPropertiesDialogDeprecated::setToAllowedList()
{
    ui->set_as_query_rb->click();
    ui->edit_query_tab_wdg->setCurrentWidget(ui->general_tab);
    setQueryTypeComboValue(QueryType::ALLOWED_LIST);
}

static void updateAtomicNumber(sketcherAtom& atom, std::string symbol)
{
    auto atomic_number = symbol_to_atomic_number(symbol);
    atom.setAtomType(atomic_number);
}

static std::set<unsigned int> parseElementList(const std::string& symbol_list)
{
    std::stringstream ss(symbol_list);
    std::set<unsigned int> atomic_nums;
    std::string symbol;
    while (ss.good()) {
        getline(ss, symbol, ',');
        auto atomic_number = symbol_to_atomic_number(symbol);
        atomic_nums.insert(atomic_number);
    }
    return atomic_nums;
}

void EditAtomPropertiesDialogDeprecated::writeAtomInfo()
{
    updateAtomicNumber(m_atom, ui->element_le->text().toStdString());
    ui->atom_common_props_wdg->writeAtomInfo(m_atom);

    auto rdk_atom = m_atom.getRDKAtom();
    if (rdk_atom != nullptr && rdk_atom->hasQuery()) {
        rdk_atom->setQuery(nullptr);
    }
}

void EditAtomPropertiesDialogDeprecated::writeQueryInfo()
{
    switch (getQueryTypeComboValue()) {
        case QueryType::SPECIFIC_ELEMENT: {
            updateAtomicNumber(m_atom,
                               ui->specific_element_le->text().toStdString());
            break;
        }
        case QueryType::ALLOWED_LIST: {
            auto atomic_nums =
                parseElementList(ui->element_list_le->text().toStdString());
            if (!atomic_nums.empty()) {
                m_atom.setAsAllowedAtomListQuery(atomic_nums);
            }
            break;
        }
        case QueryType::NOT_ALLOWED_LIST: {
            auto atomic_nums =
                parseElementList(ui->element_list_le->text().toStdString());
            if (!atomic_nums.empty()) {
                m_atom.setAsDisallowedAtomListQuery(atomic_nums);
            }
            break;
        }
        case QueryType::WILDCARD: {
            auto atom_query = getWildcardComboValue();
            m_atom.setAtomType(get_atomtype_for_atomquery(atom_query));
            break;
        }
        case QueryType::RGROUP: {
            unsigned int rgroup_number =
                static_cast<unsigned int>(ui->rgroup_sb->value());
            m_atom.setAsRGroup(rgroup_number);
            break;
        }
    }

    // A literal SMARTS query overrides the other options
    auto smarts_query = ui->smarts_query_le->text();
    if (smarts_query.isEmpty()) {

        m_atom.unsetStringProperty(ATOM_SMARTS_KEY);

        auto total_h = ui->total_H_le->text();
        if (total_h.isEmpty()) {
            m_atom.unsetIntProperty(TOTAL_H_COUNT_KEY);
        } else {
            auto h_count = total_h.toUInt();
            m_atom.setIntProperty(TOTAL_H_COUNT_KEY, h_count);
        }

        auto num_connections = ui->num_connections_le->text();
        if (num_connections.isEmpty()) {
            m_atom.unsetIntProperty(NUMBER_OF_CONNECTIONS_KEY);
        } else {
            auto connections = num_connections.toUInt();
            m_atom.setIntProperty(NUMBER_OF_CONNECTIONS_KEY, connections);
        }

        auto aromaticity = ui->aromaticity_combo->currentText();
        if (aromaticity == AROMATIC_COMBO_KEY) {
            m_atom.setIntProperty(AROMATICITY_KEY, sketcherAtom::AROMATIC_ATOM);
        } else if (aromaticity == ALIPHATIC_COMBO_KEY) {
            m_atom.setIntProperty(AROMATICITY_KEY,
                                  sketcherAtom::ALIPHATIC_ATOM);
        } else {
            m_atom.unsetIntProperty(AROMATICITY_KEY);
        }

        auto ring_count = ui->ring_count_combo->currentText();
        if (ring_count == NONZERO_COMBO_KEY) {
            m_atom.setIntProperty(MINIMUM_NUMBER_OF_RINGS_KEY, 1);
            m_atom.setIntProperty(MAXIMUM_NUMBER_OF_RINGS_KEY,
                                  MAXIMUM_RING_NUMBER_LIMIT);
        } else if (ring_count == EXACT_COMBO_KEY) {
            auto num_rings = ui->ring_count_sb->value();
            m_atom.setIntProperty(MINIMUM_NUMBER_OF_RINGS_KEY, num_rings);
            m_atom.setIntProperty(MAXIMUM_NUMBER_OF_RINGS_KEY, num_rings);
        } else {
            m_atom.unsetIntProperty(MINIMUM_NUMBER_OF_RINGS_KEY);
            m_atom.unsetIntProperty(MAXIMUM_NUMBER_OF_RINGS_KEY);
        }

        auto ring_bond_count = ui->ring_bond_count_combo->currentText();
        if (ring_bond_count == EXACT_COMBO_KEY) {
            auto num_ring_bonds = ui->ring_bond_count_sb->value();
            m_atom.setIntProperty(RING_BOND_COUNT_KEY, num_ring_bonds);
        } else {
            m_atom.unsetIntProperty(RING_BOND_COUNT_KEY);
        }

        auto smallest_ring_size = ui->smallest_ring_size_le->text();
        if (smallest_ring_size.isEmpty()) {
            m_atom.unsetIntProperty(MINIMUM_SIZE_OF_RINGS_KEY);
            m_atom.unsetIntProperty(MAXIMUM_SIZE_OF_RINGS_KEY);
        } else {
            auto ring_size = smallest_ring_size.toUInt();
            m_atom.setIntProperty(MINIMUM_SIZE_OF_RINGS_KEY, ring_size);
            m_atom.setIntProperty(MAXIMUM_SIZE_OF_RINGS_KEY,
                                  MAXIMUM_RING_SIZE_LIMIT);
        }
    } else {

        // Make sure the SMARTS is parseable, or do nothing
        auto smarts = smarts_query.toStdString();
        std::unique_ptr<RDKit::Atom> atom{RDKit::SmartsToAtom(smarts)};
        if (atom != nullptr && atom->hasQuery()) {
            parse_rdk_atom_queries(&m_atom, atom.get());

            // if we parsed a SMARTS, skip writing the info from the other
            // tab, since the SMARTS might have overwritten some of the info
            return;
        }
    }

    ui->query_common_props_wdg->writeAtomInfo(m_atom);
}

void EditAtomPropertiesDialogDeprecated::updateOKButtonEnabled()
{
    bool enable = false;
    if (ui->set_as_query_rb->isChecked()) {
        // Atom object will be set as a query. Some queries require further
        // input.
        switch (getQueryTypeComboValue()) {
            case QueryType::SPECIFIC_ELEMENT: {
                enable = !ui->specific_element_le->text().isEmpty();
                break;
            }
            case QueryType::ALLOWED_LIST:
            case QueryType::NOT_ALLOWED_LIST: {
                enable = !ui->element_list_le->text().isEmpty();
                break;
            }
            case QueryType::RGROUP:
            case QueryType::WILDCARD: {
                enable = true;
                break;
            }
        }
    } else {
        // Atom object will be set as proper atom
        enable = !ui->element_le->text().isEmpty();
    }
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(enable);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/edit_atom_properties_deprecated.moc"
