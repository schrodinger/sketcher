#include "schrodinger/sketcher/model/sketcher_model.h"

#include <stdexcept>
#include <string>

#include <QMetaType>
#include <QPointF>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"

using MonomericNucleotide = std::tuple<QString, QString, QString>;
Q_DECLARE_METATYPE(MonomericNucleotide);

namespace schrodinger
{
namespace sketcher
{

QString std_nucleobase_to_qstring(StdNucleobase base, QString U_or_T)
{
    switch (base) {
        case StdNucleobase::A:
            return "A";
        case StdNucleobase::U_OR_T:
            return U_or_T;
        case StdNucleobase::G:
            return "G";
        case StdNucleobase::C:
            return "C";
        case StdNucleobase::N:
            return "N";
        default:
            return "";
    }
}

std::vector<ModelKey> get_model_keys()
{
    return {
        ModelKey::NEW_STRUCTURES_REPLACE_CONTENT,
        ModelKey::SHOW_LID_LEGEND,
        ModelKey::ALLOW_MULTIPLE_RXNS,
        ModelKey::SHOW_VALENCE_ERRORS,
        ModelKey::COLOR_HETEROATOMS,
        ModelKey::SHOW_STEREO_LABELS,
        ModelKey::USE_IMPLICIT_HYDROGENS,
        ModelKey::SELECTION_TOOL,
        ModelKey::DRAW_TOOL,
        ModelKey::ATOM_TOOL,
        ModelKey::BOND_TOOL,
        ModelKey::CHARGE_TOOL,
        ModelKey::RING_TOOL,
        ModelKey::ENUMERATION_TOOL,
        ModelKey::ELEMENT,
        ModelKey::ATOM_QUERY,
        ModelKey::RGROUP_NUMBER,
        ModelKey::RESIDUE_TYPE,
        ModelKey::MONOMER_TOOL_TYPE,
        ModelKey::AMINO_ACID_TOOL,
        ModelKey::NUCLEIC_ACID_TOOL,
        ModelKey::RNA_NUCLEOBASE,
        ModelKey::DNA_NUCLEOBASE,
        ModelKey::CUSTOM_NUCLEOTIDE,
        ModelKey::INTERFACE_TYPE,
        ModelKey::TOOL_SET,
        ModelKey::MOLECULE_TYPE,
    };
}

SketcherModel::SketcherModel(QObject* parent) : QObject(parent)
{
    // Initialize model to default state
    reset();
    connect(this, &SketcherModel::selectionChanged, this,
            &SketcherModel::onSelectionChanged);
    connect(this, &SketcherModel::interactiveItemsChanged, this,
            &SketcherModel::onInteractiveItemsChanged);
}

QVariant SketcherModel::getValue(ModelKey key) const
{
    return m_model_map.at(key);
}

SelectionTool SketcherModel::getSelectionTool() const
{
    return m_model_map.at(ModelKey::SELECTION_TOOL).value<SelectionTool>();
}

DrawTool SketcherModel::getDrawTool() const
{
    return m_model_map.at(ModelKey::DRAW_TOOL).value<DrawTool>();
}

AtomTool SketcherModel::getAtomTool() const
{
    return m_model_map.at(ModelKey::ATOM_TOOL).value<AtomTool>();
}

BondTool SketcherModel::getBondTool() const
{
    return m_model_map.at(ModelKey::BOND_TOOL).value<BondTool>();
}

ChargeTool SketcherModel::getChargeTool() const
{
    return m_model_map.at(ModelKey::CHARGE_TOOL).value<ChargeTool>();
}

RingTool SketcherModel::getRingTool() const
{
    return m_model_map.at(ModelKey::RING_TOOL).value<RingTool>();
}

EnumerationTool SketcherModel::getEnumerationTool() const
{
    return m_model_map.at(ModelKey::ENUMERATION_TOOL).value<EnumerationTool>();
}

Element SketcherModel::getElement() const
{
    return m_model_map.at(ModelKey::ELEMENT).value<Element>();
}

AtomQuery SketcherModel::getAtomQuery() const
{
    return m_model_map.at(ModelKey::ATOM_QUERY).value<AtomQuery>();
}

MonomerToolType SketcherModel::getMonomerToolType() const
{
    return m_model_map.at(ModelKey::MONOMER_TOOL_TYPE).value<MonomerToolType>();
}

AminoAcidTool SketcherModel::getAminoAcidTool() const
{
    return m_model_map.at(ModelKey::AMINO_ACID_TOOL).value<AminoAcidTool>();
}

NucleicAcidTool SketcherModel::getNucleicAcidTool() const
{
    return m_model_map.at(ModelKey::NUCLEIC_ACID_TOOL).value<NucleicAcidTool>();
}

StdNucleobase SketcherModel::getRNANucleobase() const
{
    return m_model_map.at(ModelKey::RNA_NUCLEOBASE).value<StdNucleobase>();
}

StdNucleobase SketcherModel::getDNANucleobase() const
{
    return m_model_map.at(ModelKey::DNA_NUCLEOBASE).value<StdNucleobase>();
}

std::tuple<QString, QString, QString> SketcherModel::getCustomNucleotide() const
{
    return m_model_map.at(ModelKey::CUSTOM_NUCLEOTIDE)
        .value<MonomericNucleotide>();
}

std::optional<std::tuple<QString, QString, QString>>
SketcherModel::getNucleotide() const
{
    if (getToolSet() != ToolSet::MONOMERIC ||
        getMonomerToolType() != MonomerToolType::NUCLEIC_ACID) {
        return std::nullopt;
    } else if (getNucleicAcidTool() == NucleicAcidTool::RNA_NUCLEOTIDE) {
        auto base = std_nucleobase_to_qstring(getRNANucleobase(), "U");
        return {{"R", base, "P"}};
    } else if (getNucleicAcidTool() == NucleicAcidTool::DNA_NUCLEOTIDE) {
        auto base = std_nucleobase_to_qstring(getDNANucleobase(), "T");
        return {{"dR", base, "P"}};
    } else if (getNucleicAcidTool() == NucleicAcidTool::CUSTOM_NUCLEOTIDE) {
        return getCustomNucleotide();
    } else {
        return std::nullopt;
    }
}

InterfaceTypeType SketcherModel::getInterfaceType() const
{
    return m_model_map.at(ModelKey::INTERFACE_TYPE).value<InterfaceTypeType>();
}

ToolSet SketcherModel::getToolSet() const
{
    return m_model_map.at(ModelKey::TOOL_SET).value<ToolSet>();
}

MoleculeType SketcherModel::getMoleculeType() const
{
    return m_model_map.at(ModelKey::MOLECULE_TYPE).value<MoleculeType>();
}

bool SketcherModel::getValueBool(ModelKey key) const
{
    return m_model_map.at(key).value<bool>();
}

int SketcherModel::getValueInt(ModelKey key) const
{
    return m_model_map.at(key).value<int>();
}

QString SketcherModel::getValueString(ModelKey key) const
{
    return m_model_map.at(key).value<QString>();
}

void SketcherModel::reset()
{
    // Reset all model values to defaults
    m_model_map = {
        {ModelKey::NEW_STRUCTURES_REPLACE_CONTENT, true},
        {ModelKey::SHOW_LID_LEGEND, false},
        {ModelKey::ALLOW_MULTIPLE_RXNS, false},
        {ModelKey::SHOW_VALENCE_ERRORS, true},
        {ModelKey::COLOR_HETEROATOMS, true},
        {ModelKey::SHOW_STEREO_LABELS, true},
        {ModelKey::USE_IMPLICIT_HYDROGENS, false},
        {ModelKey::SELECTION_TOOL,
         QVariant::fromValue(SelectionTool::RECTANGLE)},
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ATOM)},
        {ModelKey::ATOM_TOOL, QVariant::fromValue(AtomTool::ELEMENT)},
        {ModelKey::BOND_TOOL, QVariant::fromValue(BondTool::SINGLE)},
        {ModelKey::CHARGE_TOOL, QVariant::fromValue(ChargeTool::DECREASE)},
        {ModelKey::RING_TOOL, QVariant::fromValue(RingTool::CYCLOPROPANE)},
        {ModelKey::ENUMERATION_TOOL,
         QVariant::fromValue(EnumerationTool::NEW_RGROUP)},
        {ModelKey::ELEMENT, QVariant::fromValue(Element::C)},
        {ModelKey::ATOM_QUERY, QVariant::fromValue(AtomQuery::A)},
        {ModelKey::RGROUP_NUMBER, 1u},
        {ModelKey::RESIDUE_TYPE, QString("")},
        {ModelKey::MONOMER_TOOL_TYPE,
         QVariant::fromValue(MonomerToolType::AMINO_ACID)},
        {ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(AminoAcidTool::ALA)},
        {ModelKey::NUCLEIC_ACID_TOOL,
         QVariant::fromValue(NucleicAcidTool::RNA_NUCLEOTIDE)},
        {ModelKey::RNA_NUCLEOBASE, QVariant::fromValue(StdNucleobase::A)},
        {ModelKey::DNA_NUCLEOBASE, QVariant::fromValue(StdNucleobase::A)},
        {ModelKey::CUSTOM_NUCLEOTIDE,
         QVariant::fromValue(MonomericNucleotide("R", "N", "P"))},
        {ModelKey::INTERFACE_TYPE, InterfaceType::ATOMISTIC},
        {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::ATOMISTIC)},
        {ModelKey::MOLECULE_TYPE, QVariant::fromValue(MoleculeType::EMPTY)},
    };
    // Notify listeners that all values have changed
    std::unordered_set<ModelKey> all_keys;
    all_keys.reserve(m_model_map.size());
    for (const auto& [key, _] : m_model_map) {
        all_keys.insert(key);
    }
    emit valuesChanged(all_keys);
}

void SketcherModel::setValues(
    const std::unordered_map<ModelKey, QVariant>& key_value_map)
{
    // Assign model values individually
    std::unordered_set<ModelKey> changed_keys;
    for (const auto& [key, value] : key_value_map) {
        auto current_value = getValue(key);
        // Outright forbid setting a value of a different type
        if (current_value.typeId() != value.typeId()) {
            throw std::runtime_error(std::string("ModelKey must be type ") +
                                     current_value.typeName());
        }
        if (current_value != value) {
            m_model_map[key] = value;
            changed_keys.insert(key);
        }
    }
    updateAtomDisplaySettings(changed_keys);

    // Finally, emit necessary signals all at once
    for (const auto& [key, value] : key_value_map) {
        emit valuePinged(key, value);
    }
    if (changed_keys.size() > 0) {
        emit valuesChanged(changed_keys);
    }
}

/**
 * @return true if the given color scheme has separate colors for heteroatoms,
 * or false if the color scheme uses the same colors for all atom types
 */
static bool does_scheme_have_heteroatoms_colored(const ColorScheme scheme)
{
    return scheme != ColorScheme::BLACK_WHITE &&
           scheme != ColorScheme::WHITE_BLACK;
}

/**
 * For a given "color heteratoms" setting, determine the "equivalent" color
 * scheme to the one that's passed in. If the passed-in color scheme matches the
 * color_heteroatoms settings, it will be returned unchanged. Otherwise, this
 * function will return a color scheme with the same background color (white or
 * black) as the passed-in scheme.
 */
static ColorScheme
get_scheme_for_color_heteroatoms_value(const ColorScheme scheme,
                                       const bool color_heteroatoms)
{
    if (color_heteroatoms == does_scheme_have_heteroatoms_colored(scheme)) {
        // the specified scheme already matches the passed color_heteroatoms
        // setting
        return scheme;
    }
    switch (scheme) {
        case ColorScheme::DEFAULT:
        case ColorScheme::AVALON:
        case ColorScheme::CDK:
            return ColorScheme::BLACK_WHITE;
        case ColorScheme::DARK_MODE:
            return ColorScheme::WHITE_BLACK;
        case ColorScheme::BLACK_WHITE:
            return ColorScheme::DEFAULT;
        case ColorScheme::WHITE_BLACK:
            return ColorScheme::DARK_MODE;
        default:
            return ColorScheme::DEFAULT;
    };
}

void SketcherModel::updateAtomDisplaySettings(
    const std::unordered_set<ModelKey>& keys)
{
    bool emit_display_settings_changed = false;
    for (auto key : keys) {
        switch (key) {
            case ModelKey::COLOR_HETEROATOMS: {
                bool color_heteroatoms = getValueBool(key);
                auto new_scheme = get_scheme_for_color_heteroatoms_value(
                    m_color_scheme, color_heteroatoms);
                setColorScheme(new_scheme);
                break;
            }
            case ModelKey::SHOW_STEREO_LABELS: {
                auto stereo_labels_shown = getValueBool(key);
                m_atom_display_settings.m_stereo_labels_visibility =
                    (stereo_labels_shown ? StereoLabels::ALL
                                         : StereoLabels::NONE);
                m_bond_display_settings.m_stereo_labels_shown =
                    stereo_labels_shown;
                emit_display_settings_changed = true;
                break;
            }
            case ModelKey::SHOW_VALENCE_ERRORS:
                m_atom_display_settings.m_valence_errors_shown =
                    getValueBool(key);
                emit_display_settings_changed = true;
                break;
            default:
                break;
        }
    }
    if (emit_display_settings_changed) {
        emit displaySettingsChanged();
    }
}

bool SketcherModel::hasReaction() const
{
    unsigned int reaction_count = emit reactionCountRequested();
    return reaction_count > 0;
}

void SketcherModel::onReactionArrowAdded()
{
    if (!getValueBool(ModelKey::ALLOW_MULTIPLE_RXNS) &&
        getEnumerationTool() == EnumerationTool::RXN_ARROW) {
        setValue(ModelKey::ENUMERATION_TOOL, EnumerationTool::RXN_PLUS);
    }
}

bool SketcherModel::sceneIsEmpty() const
{
    return getInteractiveItems().isEmpty();
}

bool SketcherModel::hasActiveSelection() const
{
    return !getSelection().isEmpty();
}

namespace
{
template <typename T> bool contains_item(const SketcherModel& model)
{
    auto contains = [](QGraphicsItem* item) {
        return dynamic_cast<T*>(item) != nullptr;
    };
    auto selection = model.getSelection();
    return std::any_of(selection.begin(), selection.end(), contains);
}
} // namespace

bool SketcherModel::hasAtomSelection() const
{
    return contains_item<AbstractAtomOrMonomerItem>(*this);
}

bool SketcherModel::hasBondSelection() const
{
    return contains_item<BondItem>(*this);
}

bool SketcherModel::hasNonMolecularObjectSelection() const
{
    return contains_item<NonMolecularItem>(*this);
}

QList<QGraphicsItem*> SketcherModel::getInteractiveItems() const
{
    return emit interactiveItemsRequested();
}

QList<QGraphicsItem*> SketcherModel::getSelection() const
{
    return emit selectionRequested();
}

QList<QGraphicsItem*> SketcherModel::getContextMenuObjects() const
{
    return emit contextMenuObjectsRequested();
}

std::pair<bool, bool> SketcherModel::getUndoStackData() const
{
    bool can_undo = emit undoStackCanUndoRequested();
    bool can_redo = emit undoStackCanRedoRequested();
    return {can_undo, can_redo};
}

std::set<unsigned int> SketcherModel::getRGroupNumbers() const
{
    std::set<unsigned int> rgroup_numbers;
    emit rGroupNumbersRequested(rgroup_numbers);
    return rgroup_numbers;
}

unsigned int SketcherModel::getNextRGroupNumber() const
{
    auto rgroup_numbers = getRGroupNumbers();
    unsigned int idx = 1;
    for (idx = 1; idx <= rgroup_numbers.size(); ++idx) {
        if (rgroup_numbers.count(idx) == 0) {
            break;
        }
    }
    return idx;
}

void SketcherModel::onSelectionChanged()
{
    // If there's a selection present, make sure that either the select tool or
    // the move-rotate tool is equipped. If not, switch to the last-used select
    // tool.
    auto draw_tool = getDrawTool();
    if (hasActiveSelection() && draw_tool != DrawTool::SELECT &&
        draw_tool != DrawTool::MOVE_ROTATE) {
        setValue(ModelKey::DRAW_TOOL, DrawTool::SELECT);
    }
}

bool SketcherModel::allItemsSelected() const
{
    auto selection = getSelection();
    auto contents = getInteractiveItems();
    return std::unordered_set<QGraphicsItem*>{selection.begin(),
                                              selection.end()} ==
           std::unordered_set<QGraphicsItem*>{contents.begin(), contents.end()};
}

void SketcherModel::onInteractiveItemsChanged()
{
    if (!sceneIsEmpty() || m_allow_select_tool_when_scene_empty) {
        return;
    }
    // If the scene is empty, select, move-rotate and delete tools should be
    // disabled. Switch to a default (draw C) tool instead
    auto draw_tool = getDrawTool();

    if ((draw_tool == DrawTool::SELECT || draw_tool == DrawTool::MOVE_ROTATE ||
         draw_tool == DrawTool::ERASE)) {
        std::unordered_map<ModelKey, QVariant> kv_pairs = {
            {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ATOM)},
            {ModelKey::ATOM_TOOL, QVariant::fromValue(AtomTool::ELEMENT)},
            {ModelKey::ELEMENT, QVariant::fromValue(Element::C)},
        };
        setValues(kv_pairs);
    }
}

void SketcherModel::setAtomDisplaySettings(
    const AtomDisplaySettings& atom_display_settings)
{
    m_atom_display_settings = atom_display_settings;
    emit displaySettingsChanged();
}

void SketcherModel::setBondDisplaySettings(
    const BondDisplaySettings& bond_display_settings)
{
    m_bond_display_settings = bond_display_settings;
    emit displaySettingsChanged();
}
void SketcherModel::setBondWidthScale(qreal scale)
{
    m_bond_display_settings.setScale(scale);
    emit displaySettingsChanged();
}

void SketcherModel::setColorScheme(const ColorScheme scheme)
{
    m_color_scheme = scheme;
    setBackgroundColor(hasDarkColorScheme() ? DARK_BACKGROUND_COLOR
                                            : LIGHT_BACKGROUND_COLOR);
    m_atom_display_settings.setColorScheme(scheme);
    m_bond_display_settings.setColorScheme(scheme);
    // if this color scheme changes the COLOR_HETEROATOMS setting, update the
    // setting and notify any listeners
    auto new_color_heteroatoms = does_scheme_have_heteroatoms_colored(scheme);
    if (new_color_heteroatoms != m_model_map[ModelKey::COLOR_HETEROATOMS]) {
        m_model_map[ModelKey::COLOR_HETEROATOMS] = new_color_heteroatoms;
        emit valuePinged(ModelKey::COLOR_HETEROATOMS, new_color_heteroatoms);
        emit valuesChanged({ModelKey::COLOR_HETEROATOMS});
    }
    emit displaySettingsChanged();
}

ColorScheme SketcherModel::getColorScheme() const
{
    return m_color_scheme;
}

std::pair<ColorScheme, ColorScheme> SketcherModel::getColorSchemes() const
{
    return {
        get_scheme_for_color_heteroatoms_value(m_color_scheme, true),
        get_scheme_for_color_heteroatoms_value(m_color_scheme, false),
    };
}

void SketcherModel::setBackgroundColor(QColor color)
{
    m_background_color = color;
    emit backgroundColorChanged(color);
}

bool SketcherModel::hasDarkColorScheme() const
{
    return DARK_MODE_COLOR_SCHEMES.contains(m_color_scheme);
}

void SketcherModel::setFontSize(int size)
{
    m_font_size = size;
    emit displaySettingsChanged();
}

int SketcherModel::getFontSize() const
{
    // Apply debug mode scaling if atom indices are shown
    constexpr double DEBUG_MODE_FONT_SCALING = 0.7;
    if (m_atom_display_settings.m_show_atom_indices) {
        return static_cast<int>(m_font_size * DEBUG_MODE_FONT_SCALING);
    }
    return m_font_size;
}

void SketcherModel::loadRenderOptions(const RenderOptions& opts)
{
    m_font_size = opts.font_size;
    AtomDisplaySettings atom_display_settings(m_atom_display_settings);
    atom_display_settings.m_carbon_labels = opts.carbon_labels;
    atom_display_settings.m_show_symbol_for_H_isotopes =
        opts.show_symbol_for_H_isotopes;
    atom_display_settings.m_explicit_abs_labels_shown =
        opts.show_absolute_stereo_groups;
    atom_display_settings.m_show_simplified_stereo_annotation =
        opts.show_simplified_stereo_annotation;
    atom_display_settings.m_stereo_labels_visibility =
        opts.show_stereo_annotations;
    atom_display_settings.setSquigglePenScale(opts.bond_width_scale);
    atom_display_settings.setColorScheme(opts.color_scheme);
    setAtomDisplaySettings(atom_display_settings);

    BondDisplaySettings bond_display_settings(m_bond_display_settings);
    bond_display_settings.setScale(opts.bond_width_scale);
    bond_display_settings.m_color =
        atom_display_settings.getAtomColor(static_cast<int>(Element::C));
    bond_display_settings.m_stereo_labels_shown =
        opts.show_stereo_annotations != StereoLabels::NONE;
    setBondDisplaySettings(bond_display_settings);
}

void SketcherModel::setSelectToolAllowedWhenSceneEmpty(const bool allowed)
{
    m_allow_select_tool_when_scene_empty = allowed;
}

void SketcherModel::setSelectOnlyModeActive(const bool select_only_mode)
{
    m_select_only_mode_active = select_only_mode;
}

bool SketcherModel::isSelectOnlyModeActive() const
{
    return m_select_only_mode_active;
}

const AtomDisplaySettings* SketcherModel::getAtomDisplaySettingsPtr() const
{
    return &m_atom_display_settings;
}

const BondDisplaySettings* SketcherModel::getBondDisplaySettingsPtr() const
{
    return &m_bond_display_settings;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/model/sketcher_model.moc"
