/**
 * Qt-free filtering, searching and sorting logic for the Monomer Database
 * popup. The Qt table model and popup widgets are thin adapters over these
 * classes and functions.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Bit flags for the polymer types that the Monomer Database popup can filter
 * on. NUCLEIC_ACID corresponds to rdkit_extensions::ChainType::RNA (which also
 * covers DNA monomers) and CHEM is shown to the user as "Structure".
 */
enum class MonomerDbPolymerType : unsigned {
    NONE = 0,
    PEPTIDE = 1 << 0,
    NUCLEIC_ACID = 1 << 1,
    CHEM = 1 << 2,
    ALL = PEPTIDE | NUCLEIC_ACID | CHEM,
};

constexpr MonomerDbPolymerType operator|(MonomerDbPolymerType lhs,
                                         MonomerDbPolymerType rhs)
{
    return static_cast<MonomerDbPolymerType>(static_cast<unsigned>(lhs) |
                                             static_cast<unsigned>(rhs));
}

constexpr MonomerDbPolymerType operator&(MonomerDbPolymerType lhs,
                                         MonomerDbPolymerType rhs)
{
    return static_cast<MonomerDbPolymerType>(static_cast<unsigned>(lhs) &
                                             static_cast<unsigned>(rhs));
}

/**
 * @return whether all of the flags in `flag` are set in `polymer_types`. Note
 * that this returns false if `flag` is NONE.
 */
constexpr bool has_flag(MonomerDbPolymerType polymer_types,
                        MonomerDbPolymerType flag)
{
    return flag != MonomerDbPolymerType::NONE && (polymer_types & flag) == flag;
}

/**
 * @return the number of flags set in `polymer_types`
 */
constexpr int count_flags(MonomerDbPolymerType polymer_types)
{
    int count = 0;
    for (auto flag :
         {MonomerDbPolymerType::PEPTIDE, MonomerDbPolymerType::NUCLEIC_ACID,
          MonomerDbPolymerType::CHEM}) {
        if (has_flag(polymer_types, flag)) {
            ++count;
        }
    }
    return count;
}

/**
 * Values for the Monomer type filter, which only applies to peptides
 */
enum class MonomerTypeFilter {
    ANY,
    BACKBONE,
    TERMINAL,
};

/**
 * Values for the Natural analog filter
 */
enum class NaturalAnalogFilter {
    ANY_OR_NONE, ///< no analog filtering
    ANY,         ///< monomers that have a natural analog
    NONE,        ///< monomers without a natural analog (see
                 ///< has_no_natural_analog)
    SPECIFIC,    ///< monomers with a specific natural analog
};

/**
 * The columns of the Monomer Database table, in display order
 */
enum class MonomerDbColumn {
    SYMBOL,
    NAME,
    POLYMER_TYPE,
    NATURAL_ANALOG,
};

/**
 * The header label of each column of the Monomer Database table, in
 * MonomerDbColumn order. Each label is paired with its column so that the
 * static_assert below can check the order: reordering MonomerDbColumn without
 * reordering this array fails to compile. A column added to MonomerDbColumn
 * must also be added here, since this array's size is the table's column
 * count.
 */
constexpr std::array<std::pair<MonomerDbColumn, std::string_view>, 4>
    MONOMER_DB_COLUMN_LABELS{{
        {MonomerDbColumn::SYMBOL, "Symbol"},
        {MonomerDbColumn::NAME, "Name"},
        {MonomerDbColumn::POLYMER_TYPE, "Type"},
        {MonomerDbColumn::NATURAL_ANALOG, "Analog"},
    }};

static_assert(
    [] {
        for (std::size_t i = 0; i < MONOMER_DB_COLUMN_LABELS.size(); ++i) {
            if (static_cast<std::size_t>(MONOMER_DB_COLUMN_LABELS[i].first) !=
                i) {
                return false;
            }
        }
        return true;
    }(),
    "MONOMER_DB_COLUMN_LABELS must list the columns in MonomerDbColumn order");

/**
 * The sort column and direction for the Monomer Database table. This is kept
 * separate from MonomerDatabaseFilter so that resetting the filter never
 * affects the sort.
 */
struct MonomerDbSort {
    MonomerDbColumn column = MonomerDbColumn::SYMBOL;
    bool ascending = true;

    bool operator==(const MonomerDbSort&) const = default;
};

/**
 * The filter settings for the Monomer Database popup.
 *
 * The Monomer type and Natural analog settings only apply for certain polymer
 * types (see isMonomerTypeApplicable() and isNaturalAnalogApplicable()). While
 * a setting is inapplicable, its effective value (as returned by the getters
 * and as used by matches(), the filter count and the description) is the most
 * permissive value, but the user's choice is remembered and restored once the
 * setting applies again. The exceptions are:
 *  - the remembered natural analog is forgotten when the setting applies again
 *    for a different polymer type (Peptide vs. Nucleic acid) than the one it
 *    was chosen for
 *  - the remembered monomer type is forgotten once the natural analog applies
 *    to Nucleic acid
 *
 * A default-constructed filter is the most permissive filter: all polymer
 * types, monomer type ANY and natural analog ANY_OR_NONE. "Filtering" (i.e.
 * the filter count and the is*Filtered() methods) is relative to that filter.
 */
class SKETCHER_API MonomerDatabaseFilter
{
  public:
    MonomerDatabaseFilter() = default;

    /**
     * Create a filter for the given polymer types, with monomer type ANY and
     * natural analog ANY_OR_NONE.
     */
    explicit MonomerDatabaseFilter(MonomerDbPolymerType polymer_types);

    /**
     * Set the polymer types to show.
     *
     * Note that this may clear a remembered setting: the remembered natural
     * analog and monomer type will be forgotten when the natural analog setting
     * becomes applicable for a different polymer type than the one the value
     * was chosen for. For example, switching filtering from peptides to both
     * peptides and nucleic acids will remember the natural analog and monomer
     * settings from peptide (so that switching the polymer type filtering back
     * to peptide will restore them). However, switching from peptides to both
     * peptides and nucleic acids and then switching again to only nucleic acids
     * will clear the remembered natural analog and monomer settings.
     */
    void setPolymerTypes(MonomerDbPolymerType polymer_types);

    /**
     * Set the monomer type. The value is remembered even if the monomer type
     * isn't currently applicable, but it only affects the filtering while it
     * is applicable.
     */
    void setMonomerType(MonomerTypeFilter monomer_type);

    /**
     * Set the natural analog, and record the currently applicable polymer type
     * (Peptide or Nucleic acid) as the type that it was chosen for.
     * @param analog The natural analog setting
     * @param specific The natural analog symbol to match. Only used if
     * `analog` is SPECIFIC.
     * @pre isNaturalAnalogApplicable()
     */
    void setNaturalAnalog(NaturalAnalogFilter analog,
                          std::string specific = {});

    /**
     * @return the polymer types to show
     */
    MonomerDbPolymerType getPolymerTypes() const;

    /**
     * @return the effective monomer type, which is ANY while the monomer type
     * is inapplicable
     */
    MonomerTypeFilter getMonomerType() const;

    /**
     * @return the effective natural analog setting, which is ANY_OR_NONE while
     * the natural analog is inapplicable
     */
    NaturalAnalogFilter getNaturalAnalog() const;

    /**
     * @return the effective natural analog symbol to match, or an empty string
     * unless getNaturalAnalog() is SPECIFIC
     */
    const std::string& getSpecificNaturalAnalog() const;

    /**
     * @return whether the monomer type setting applies, which is true when
     * Peptide is shown and Nucleic acid isn't. CHEM has no effect.
     */
    bool isMonomerTypeApplicable() const;

    /**
     * @return whether the natural analog setting applies, which is true when
     * exactly one of Peptide and Nucleic acid is shown. CHEM has no effect.
     */
    bool isNaturalAnalogApplicable() const;

    /**
     * @return the specific natural analog symbols that the user can choose
     * from, which depend on the polymer type that the natural analog setting
     * applies to. Returns an empty list if the natural analog is inapplicable.
     */
    std::vector<std::string> getNaturalAnalogChoices() const;

    /**
     * @return whether the polymer types are anything other than all three
     */
    bool isPolymerTypeFiltered() const;

    /**
     * @return whether the effective monomer type is anything other than ANY
     */
    bool isMonomerTypeFiltered() const;

    /**
     * @return whether the effective natural analog is anything other than
     * ANY_OR_NONE
     */
    bool isNaturalAnalogFiltered() const;

    /**
     * @return the number of settings that are filtered, i.e. that differ from
     * a default-constructed filter
     */
    int getNumFiltersApplied() const;

    /**
     * @return whether any setting is filtered
     */
    bool isFiltering() const;

    /**
     * @return a description of the monomers that this filter shows, e.g.
     * "Peptide, backbone, Cys analogs"
     */
    std::string getDescription() const;

    /**
     * @return whether the given monomer passes this filter
     */
    bool matches(const rdkit_extensions::MonomerInfo& info) const;

    /**
     * Compares the full state, including remembered settings that are not
     * currently applicable. To check whether two filters apply the same
     * filtering, compare the getters instead.
     */
    bool operator==(const MonomerDatabaseFilter&) const = default;

  private:
    /**
     * @return the polymer type (PEPTIDE or NUCLEIC_ACID) that the natural
     * analog setting applies to, or NONE if it is inapplicable
     */
    MonomerDbPolymerType getNaturalAnalogPolymerType() const;

    MonomerDbPolymerType m_polymer_types = MonomerDbPolymerType::ALL;
    // the remembered settings, which may not be applicable
    MonomerTypeFilter m_monomer_type = MonomerTypeFilter::ANY;
    NaturalAnalogFilter m_natural_analog = NaturalAnalogFilter::ANY_OR_NONE;
    std::string m_specific_natural_analog;
    // the polymer type (PEPTIDE or NUCLEIC_ACID) that m_natural_analog was
    // chosen for, or NONE if it was never set
    MonomerDbPolymerType m_natural_analog_polymer_type =
        MonomerDbPolymerType::NONE;
};

/**
 * @return the polymer type flag (PEPTIDE, NUCLEIC_ACID or CHEM) for the given
 * monomer, or NONE if its polymer type is missing or not one of those three
 */
SKETCHER_API MonomerDbPolymerType
get_monomer_db_polymer_type(const rdkit_extensions::MonomerInfo& info);

/**
 * @return the text to display for a single polymer type flag ("Peptide",
 * "Nucleic acid" or "Structure"), or an empty string for any other value
 */
SKETCHER_API const std::string&
get_polymer_type_display_text(MonomerDbPolymerType polymer_type);

/**
 * @return whether the given monomer has no natural analog. This is true if
 * its natural analog is missing or empty, or is the unknown residue for its
 * polymer type: X for peptides or N for nucleic acids. (Note that N is
 * asparagine for peptides.) Always true for CHEM monomers.
 */
SKETCHER_API bool
has_no_natural_analog(const rdkit_extensions::MonomerInfo& info);

/**
 * @return the natural analog text to display for the given monomer, which is
 * blank if has_no_natural_analog() is true
 */
SKETCHER_API std::string
get_natural_analog_display_text(const rdkit_extensions::MonomerInfo& info);

/**
 * Build the source list for the Monomer Database table: every custom monomer,
 * plus every core monomer with no natural analog (see has_no_natural_analog)
 * unless a custom monomer has the same symbol and polymer type.
 * @param custom_monomers The custom monomers, e.g. from
 * MonomerDatabase::getCustomMonomers()
 * @param core_monomers The core monomers, e.g. from
 * MonomerDatabase::getCoreMonomers()
 * @return the custom monomers in their original order, followed by the
 * included core monomers in their original order
 */
SKETCHER_API std::vector<rdkit_extensions::MonomerInfo>
get_monomer_database_table_monomers(
    std::vector<rdkit_extensions::MonomerInfo> custom_monomers,
    const std::vector<rdkit_extensions::MonomerInfo>& core_monomers);

/**
 * Filter and sort the given monomers. The sort is case-insensitive and
 * stable, and ties are broken by symbol and then by polymer type. Missing
 * values (including blank natural analogs) sort last in ascending order and
 * first in descending order.
 * @return the indices into `all_monomers` of the matching monomers, in sorted
 * order
 */
SKETCHER_API std::vector<std::size_t> filter_and_sort_monomers(
    const std::vector<rdkit_extensions::MonomerInfo>& all_monomers,
    const MonomerDatabaseFilter& filter, const MonomerDbSort& sort);

/**
 * @param all_monomers The monomers that `rows` indexes into
 * @param rows Indices into `all_monomers`
 * @param text The search text
 * @return the subset of `rows` whose monomer's symbol or name contains `text`,
 * case-insensitively, in their original order. If `text` is empty, `rows` is
 * returned unchanged.
 */
SKETCHER_API std::vector<std::size_t> filter_by_search_text(
    const std::vector<rdkit_extensions::MonomerInfo>& all_monomers,
    const std::vector<std::size_t>& rows, std::string_view text);

} // namespace sketcher
} // namespace schrodinger
