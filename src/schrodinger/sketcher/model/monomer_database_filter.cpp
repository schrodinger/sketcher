#include "schrodinger/sketcher/model/monomer_database_filter.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <compare>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

using schrodinger::rdkit_extensions::ChainType;
using schrodinger::rdkit_extensions::MonomerInfo;

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * A matching monomer's index along with the values used to sort it. A missing
 * value is std::nullopt.
 */
struct SortKey {
    std::size_t index;
    std::optional<std::string_view> primary;
    std::optional<std::string_view> symbol;
    std::optional<std::string_view> polymer_type;
};

} // namespace

// The MONOMER_TYPE values that the Monomer type filter matches,
// case-insensitively
static const std::string BACKBONE_MONOMER_TYPE = "Backbone";
static const std::string TERMINAL_MONOMER_TYPE = "Terminal";

static const std::string PEPTIDE_DISPLAY_TEXT = "Peptide";
static const std::string NUCLEIC_ACID_DISPLAY_TEXT = "Nucleic acid";
static const std::string STRUCTURE_DISPLAY_TEXT = "Structure";

static constexpr std::array<MonomerDbPolymerType, 3> ALL_POLYMER_TYPE_FLAGS{
    MonomerDbPolymerType::PEPTIDE, MonomerDbPolymerType::NUCLEIC_ACID,
    MonomerDbPolymerType::CHEM};

static char to_lower(char c)
{
    return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
}

static bool chars_equal_case_insensitive(char a, char b)
{
    return to_lower(a) == to_lower(b);
}

static bool equals_case_insensitive(std::string_view a, std::string_view b)
{
    return std::ranges::equal(a, b, chars_equal_case_insensitive);
}

static bool contains_case_insensitive(std::string_view haystack,
                                      std::string_view needle)
{
    return !std::ranges::search(haystack, needle, chars_equal_case_insensitive)
                .empty();
}

static std::strong_ordering compare_case_insensitive(std::string_view a,
                                                     std::string_view b)
{
    return std::lexicographical_compare_three_way(
        a.begin(), a.end(), b.begin(), b.end(),
        [](char x, char y) { return to_lower(x) <=> to_lower(y); });
}

/**
 * Compare two sort values case-insensitively, with missing values after all
 * other values
 */
static std::strong_ordering
compare_sort_values(const std::optional<std::string_view>& a,
                    const std::optional<std::string_view>& b)
{
    if (!a.has_value() || !b.has_value()) {
        // true (i.e. missing) sorts after false
        return b.has_value() <=> a.has_value();
    }
    return compare_case_insensitive(*a, *b);
}

/**
 * @return the given value as a sort value, where a missing or empty value is
 * std::nullopt
 */
static std::optional<std::string_view>
to_sort_value(const std::optional<std::string>& value)
{
    if (!value.has_value() || value->empty()) {
        return std::nullopt;
    }
    return *value;
}

static std::optional<std::string_view> to_sort_value(const std::string& value)
{
    if (value.empty()) {
        return std::nullopt;
    }
    return value;
}

/**
 * @return the symbol of the unknown residue for the given polymer type (X for
 * peptides, N for nucleic acids), or nullptr for any other polymer type
 */
static const std::string*
get_unknown_residue_symbol(MonomerDbPolymerType polymer_type)
{
    switch (polymer_type) {
        case MonomerDbPolymerType::PEPTIDE:
            return &AMINO_ACID_TOOL_TO_RES_NAME.at(AminoAcidTool::UNK);
        case MonomerDbPolymerType::NUCLEIC_ACID:
            return &NUCLEIC_ACID_TOOL_TO_RES_NAME.at(NucleicAcidTool::N);
        default:
            return nullptr;
    }
}

/**
 * Determine whether the given monomer has no natural analog, given its
 * already-determined polymer type
 */
static bool has_no_natural_analog(const MonomerInfo& info,
                                  MonomerDbPolymerType polymer_type)
{
    if (polymer_type == MonomerDbPolymerType::CHEM ||
        !info.natural_analog.has_value() || info.natural_analog->empty()) {
        return true;
    }
    auto* unknown_symbol = get_unknown_residue_symbol(polymer_type);
    return unknown_symbol != nullptr && *info.natural_analog == *unknown_symbol;
}

static const std::string& get_monomer_type_text(MonomerTypeFilter monomer_type)
{
    static const std::string empty;
    switch (monomer_type) {
        case MonomerTypeFilter::BACKBONE:
            return BACKBONE_MONOMER_TYPE;
        case MonomerTypeFilter::TERMINAL:
            return TERMINAL_MONOMER_TYPE;
        default:
            return empty;
    }
}

/**
 * @return a description of the given polymer types, e.g. "Peptide or
 * structure" or "All"
 */
static std::string
get_polymer_types_description(MonomerDbPolymerType polymer_types)
{
    auto num_types = count_flags(polymer_types);
    if (num_types == static_cast<int>(ALL_POLYMER_TYPE_FLAGS.size())) {
        return "All";
    } else if (num_types == 0) {
        return "No";
    }
    std::string description;
    for (auto flag : ALL_POLYMER_TYPE_FLAGS) {
        if (!has_flag(polymer_types, flag)) {
            continue;
        }
        auto text = get_polymer_type_display_text(flag);
        if (!description.empty()) {
            description += " or ";
            text[0] = to_lower(text[0]);
        }
        description += text;
    }
    return description;
}

/**
 * @return the text used to describe a specific natural analog, which is the
 * three-letter code (e.g. "Cys") for standard amino acids and the symbol
 * otherwise (i.e. for nucleic acids)
 */
static std::string get_specific_analog_description(MonomerDbPolymerType type,
                                                   const std::string& analog)
{
    if (type == MonomerDbPolymerType::PEPTIDE) {
        auto& codes = rdkit_extensions::AMINO_ACID_THREE_LETTER_CODES;
        if (auto it = codes.find(analog); it != codes.end()) {
            // convert e.g. "CYS" to "Cys"
            std::string code = it->second;
            std::transform(code.begin() + 1, code.end(), code.begin() + 1,
                           to_lower);
            return code;
        }
    }
    return analog;
}

static SortKey make_sort_key(const std::vector<MonomerInfo>& all_monomers,
                             std::size_t index, MonomerDbColumn column)
{
    const auto& info = all_monomers[index];
    auto polymer_type = get_monomer_db_polymer_type(info);
    SortKey key{index, std::nullopt, to_sort_value(info.symbol),
                to_sort_value(get_polymer_type_display_text(polymer_type))};
    switch (column) {
        case MonomerDbColumn::SYMBOL:
            key.primary = key.symbol;
            break;
        case MonomerDbColumn::NAME:
            key.primary = to_sort_value(info.name);
            break;
        case MonomerDbColumn::POLYMER_TYPE:
            key.primary = key.polymer_type;
            break;
        case MonomerDbColumn::NATURAL_ANALOG:
            // monomers with no natural analog display a blank analog, so
            // they're sorted as missing values
            if (!has_no_natural_analog(info, polymer_type)) {
                key.primary = to_sort_value(info.natural_analog);
            }
            break;
    }
    return key;
}

MonomerDatabaseFilter::MonomerDatabaseFilter(
    MonomerDbPolymerType polymer_types) :
    m_polymer_types(polymer_types)
{
}

void MonomerDatabaseFilter::setPolymerTypes(MonomerDbPolymerType polymer_types)
{
    m_polymer_types = polymer_types;
    auto analog_type = getNaturalAnalogPolymerType();
    if (analog_type == MonomerDbPolymerType::NONE) {
        return;
    }
    if (m_natural_analog_polymer_type != MonomerDbPolymerType::NONE &&
        m_natural_analog_polymer_type != analog_type) {
        // the natural analog was chosen for the other polymer type, so forget
        // it
        m_natural_analog = NaturalAnalogFilter::ANY_OR_NONE;
        m_specific_natural_analog.clear();
        m_natural_analog_polymer_type = MonomerDbPolymerType::NONE;
    }
    if (analog_type == MonomerDbPolymerType::NUCLEIC_ACID) {
        // the user has moved on to nucleic acids, so forget the monomer type
        m_monomer_type = MonomerTypeFilter::ANY;
    }
}

void MonomerDatabaseFilter::setMonomerType(MonomerTypeFilter monomer_type)
{
    m_monomer_type = monomer_type;
}

void MonomerDatabaseFilter::setNaturalAnalog(NaturalAnalogFilter analog,
                                             std::string specific)
{
    assert(isNaturalAnalogApplicable());
    if (!isNaturalAnalogApplicable()) {
        return;
    }
    m_natural_analog = analog;
    m_specific_natural_analog =
        analog == NaturalAnalogFilter::SPECIFIC ? std::move(specific) : "";
    m_natural_analog_polymer_type = getNaturalAnalogPolymerType();
}

MonomerDbPolymerType MonomerDatabaseFilter::getPolymerTypes() const
{
    return m_polymer_types;
}

MonomerTypeFilter MonomerDatabaseFilter::getMonomerType() const
{
    return isMonomerTypeApplicable() ? m_monomer_type : MonomerTypeFilter::ANY;
}

NaturalAnalogFilter MonomerDatabaseFilter::getNaturalAnalog() const
{
    return isNaturalAnalogApplicable() ? m_natural_analog
                                       : NaturalAnalogFilter::ANY_OR_NONE;
}

const std::string& MonomerDatabaseFilter::getSpecificNaturalAnalog() const
{
    static const std::string empty;
    return getNaturalAnalog() == NaturalAnalogFilter::SPECIFIC
               ? m_specific_natural_analog
               : empty;
}

bool MonomerDatabaseFilter::isMonomerTypeApplicable() const
{
    return has_flag(m_polymer_types, MonomerDbPolymerType::PEPTIDE) &&
           !has_flag(m_polymer_types, MonomerDbPolymerType::NUCLEIC_ACID);
}

bool MonomerDatabaseFilter::isNaturalAnalogApplicable() const
{
    return getNaturalAnalogPolymerType() != MonomerDbPolymerType::NONE;
}

MonomerDbPolymerType MonomerDatabaseFilter::getNaturalAnalogPolymerType() const
{
    bool peptide = has_flag(m_polymer_types, MonomerDbPolymerType::PEPTIDE);
    bool nucleic_acid =
        has_flag(m_polymer_types, MonomerDbPolymerType::NUCLEIC_ACID);
    if (peptide == nucleic_acid) {
        return MonomerDbPolymerType::NONE;
    }
    return peptide ? MonomerDbPolymerType::PEPTIDE
                   : MonomerDbPolymerType::NUCLEIC_ACID;
}

std::vector<std::string> MonomerDatabaseFilter::getNaturalAnalogChoices() const
{
    std::vector<std::string> choices;
    switch (getNaturalAnalogPolymerType()) {
        case MonomerDbPolymerType::PEPTIDE:
            // the unknown amino acid (X) is covered by
            // NaturalAnalogFilter::NONE
            for (auto tool : AMINO_ACID_TOOL_DISPLAY_ORDER) {
                if (tool != AminoAcidTool::UNK) {
                    choices.push_back(AMINO_ACID_TOOL_TO_RES_NAME.at(tool));
                }
            }
            break;
        case MonomerDbPolymerType::NUCLEIC_ACID:
            // a sugar whose analog is deoxyribose uses R as its natural analog,
            // so dR isn't included
            for (auto tool : {NucleicAcidTool::R, NucleicAcidTool::P}) {
                choices.push_back(NUCLEIC_ACID_TOOL_TO_RES_NAME.at(tool));
            }
            for (auto tool : NA_BASE_DISPLAY_ORDER) {
                choices.push_back(NUCLEIC_ACID_TOOL_TO_RES_NAME.at(tool));
            }
            break;
        default:
            break;
    }
    return choices;
}

bool MonomerDatabaseFilter::isPolymerTypeFiltered() const
{
    return m_polymer_types != MonomerDbPolymerType::ALL;
}

bool MonomerDatabaseFilter::isMonomerTypeFiltered() const
{
    return getMonomerType() != MonomerTypeFilter::ANY;
}

bool MonomerDatabaseFilter::isNaturalAnalogFiltered() const
{
    return getNaturalAnalog() != NaturalAnalogFilter::ANY_OR_NONE;
}

int MonomerDatabaseFilter::getNumFiltersApplied() const
{
    return static_cast<int>(isPolymerTypeFiltered()) +
           static_cast<int>(isMonomerTypeFiltered()) +
           static_cast<int>(isNaturalAnalogFiltered());
}

bool MonomerDatabaseFilter::isFiltering() const
{
    return getNumFiltersApplied() > 0;
}

std::string MonomerDatabaseFilter::getDescription() const
{
    auto monomer_type = getMonomerType();
    auto analog = getNaturalAnalog();
    auto described_types = m_polymer_types;
    if (monomer_type != MonomerTypeFilter::ANY ||
        analog == NaturalAnalogFilter::ANY ||
        analog == NaturalAnalogFilter::SPECIFIC) {
        // these settings exclude all CHEM monomers, so don't mention them
        described_types =
            described_types & (MonomerDbPolymerType::PEPTIDE |
                               MonomerDbPolymerType::NUCLEIC_ACID);
    }
    auto description = get_polymer_types_description(described_types);
    if (monomer_type != MonomerTypeFilter::ANY) {
        auto monomer_type_text = get_monomer_type_text(monomer_type);
        std::ranges::transform(monomer_type_text, monomer_type_text.begin(),
                               to_lower);
        description += ", " + monomer_type_text;
    }
    switch (analog) {
        case NaturalAnalogFilter::ANY_OR_NONE:
            description += " monomers";
            break;
        case NaturalAnalogFilter::ANY:
            description += " monomers with any natural analog";
            break;
        case NaturalAnalogFilter::NONE:
            description += " monomers with no natural analog";
            break;
        case NaturalAnalogFilter::SPECIFIC:
            description +=
                ", " +
                get_specific_analog_description(getNaturalAnalogPolymerType(),
                                                m_specific_natural_analog) +
                " analogs";
            break;
    }
    return description;
}

bool MonomerDatabaseFilter::matches(const MonomerInfo& info) const
{
    auto polymer_type = get_monomer_db_polymer_type(info);
    if (!has_flag(m_polymer_types, polymer_type)) {
        return false;
    }

    if (auto monomer_type = getMonomerType();
        monomer_type != MonomerTypeFilter::ANY) {
        if (polymer_type == MonomerDbPolymerType::CHEM ||
            !info.monomer_type.has_value() ||
            !equals_case_insensitive(*info.monomer_type,
                                     get_monomer_type_text(monomer_type))) {
            return false;
        }
    }

    switch (getNaturalAnalog()) {
        case NaturalAnalogFilter::ANY_OR_NONE:
            return true;
        case NaturalAnalogFilter::ANY:
            return !has_no_natural_analog(info, polymer_type);
        case NaturalAnalogFilter::NONE:
            return has_no_natural_analog(info, polymer_type);
        case NaturalAnalogFilter::SPECIFIC:
            return !has_no_natural_analog(info, polymer_type) &&
                   *info.natural_analog == m_specific_natural_analog;
    }
    return false;
}

MonomerDbPolymerType get_monomer_db_polymer_type(const MonomerInfo& info)
{
    if (!info.polymer_type.has_value()) {
        return MonomerDbPolymerType::NONE;
    }
    ChainType chain_type;
    try {
        chain_type = rdkit_extensions::toChainType(*info.polymer_type);
    } catch (const std::invalid_argument&) {
        return MonomerDbPolymerType::NONE;
    }
    switch (chain_type) {
        case ChainType::PEPTIDE:
            return MonomerDbPolymerType::PEPTIDE;
        case ChainType::RNA:
            return MonomerDbPolymerType::NUCLEIC_ACID;
        case ChainType::CHEM:
            return MonomerDbPolymerType::CHEM;
        default:
            return MonomerDbPolymerType::NONE;
    }
}

const std::string& get_polymer_type_display_text(MonomerDbPolymerType type)
{
    static const std::string empty;
    switch (type) {
        case MonomerDbPolymerType::PEPTIDE:
            return PEPTIDE_DISPLAY_TEXT;
        case MonomerDbPolymerType::NUCLEIC_ACID:
            return NUCLEIC_ACID_DISPLAY_TEXT;
        case MonomerDbPolymerType::CHEM:
            return STRUCTURE_DISPLAY_TEXT;
        default:
            return empty;
    }
}

bool has_no_natural_analog(const MonomerInfo& info)
{
    return has_no_natural_analog(info, get_monomer_db_polymer_type(info));
}

std::string get_natural_analog_display_text(const MonomerInfo& info)
{
    if (has_no_natural_analog(info)) {
        return {};
    }
    return *info.natural_analog;
}

std::vector<MonomerInfo> get_monomer_database_table_monomers(
    std::vector<MonomerInfo> custom_monomers,
    const std::vector<MonomerInfo>& core_monomers)
{
    // (symbol, polymer type) for each custom monomer
    std::set<std::pair<std::string, std::string>> custom_ids;
    for (const auto& info : custom_monomers) {
        custom_ids.emplace(info.symbol.value_or(""),
                           info.polymer_type.value_or(""));
    }
    auto table_monomers = std::move(custom_monomers);
    for (const auto& info : core_monomers) {
        if (has_no_natural_analog(info) &&
            !custom_ids.contains(
                {info.symbol.value_or(""), info.polymer_type.value_or("")})) {
            table_monomers.push_back(info);
        }
    }
    return table_monomers;
}

std::vector<std::size_t>
filter_and_sort_monomers(const std::vector<MonomerInfo>& all_monomers,
                         const MonomerDatabaseFilter& filter,
                         const MonomerDbSort& sort)
{
    std::vector<SortKey> keys;
    for (std::size_t i = 0; i < all_monomers.size(); ++i) {
        if (filter.matches(all_monomers[i])) {
            keys.push_back(make_sort_key(all_monomers, i, sort.column));
        }
    }
    std::ranges::stable_sort(keys, [&sort](const SortKey& a, const SortKey& b) {
        if (auto cmp = compare_sort_values(a.primary, b.primary); cmp != 0) {
            return sort.ascending ? cmp < 0 : cmp > 0;
        }
        if (auto cmp = compare_sort_values(a.symbol, b.symbol); cmp != 0) {
            return cmp < 0;
        }
        return compare_sort_values(a.polymer_type, b.polymer_type) < 0;
    });
    std::vector<std::size_t> rows;
    rows.reserve(keys.size());
    for (const auto& key : keys) {
        rows.push_back(key.index);
    }
    return rows;
}

std::vector<std::size_t>
filter_by_search_text(const std::vector<MonomerInfo>& all_monomers,
                      const std::vector<std::size_t>& rows,
                      std::string_view text)
{
    if (text.empty()) {
        return rows;
    }
    auto field_contains_text = [text](const std::optional<std::string>& field) {
        return field.has_value() && contains_case_insensitive(*field, text);
    };
    std::vector<std::size_t> matching_rows;
    for (auto row : rows) {
        const auto& info = all_monomers[row];
        if (field_contains_text(info.symbol) ||
            field_contains_text(info.name)) {
            matching_rows.push_back(row);
        }
    }
    return matching_rows;
}

} // namespace sketcher
} // namespace schrodinger
