#define BOOST_TEST_MODULE test_monomer_database_filter

#include <cstddef>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/monomer_database_filter.h"

using schrodinger::rdkit_extensions::MonomerInfo;

namespace schrodinger
{
namespace sketcher
{

MAKE_ENUM_LOGGABLE(MonomerDbPolymerType)
MAKE_ENUM_LOGGABLE(MonomerTypeFilter)
MAKE_ENUM_LOGGABLE(NaturalAnalogFilter)

using PolymerType = MonomerDbPolymerType;

static const PolymerType PEPTIDE = PolymerType::PEPTIDE;
static const PolymerType NUCLEIC_ACID = PolymerType::NUCLEIC_ACID;
static const PolymerType CHEM = PolymerType::CHEM;

/**
 * Every combination of polymer type flags, from NONE to ALL
 */
static const std::vector<PolymerType> ALL_POLYMER_TYPE_COMBINATIONS{
    PolymerType::NONE,      PEPTIDE,
    NUCLEIC_ACID,           CHEM,
    PEPTIDE | NUCLEIC_ACID, PEPTIDE | CHEM,
    NUCLEIC_ACID | CHEM,    PolymerType::ALL,
};

static MonomerInfo
make_monomer(std::optional<std::string> symbol,
             std::optional<std::string> polymer_type,
             std::optional<std::string> natural_analog = std::nullopt,
             std::optional<std::string> monomer_type = std::nullopt,
             std::optional<std::string> name = std::nullopt)
{
    MonomerInfo info;
    info.symbol = std::move(symbol);
    info.polymer_type = std::move(polymer_type);
    info.natural_analog = std::move(natural_analog);
    info.monomer_type = std::move(monomer_type);
    info.name = std::move(name);
    return info;
}

static std::vector<std::string>
get_symbols(const std::vector<MonomerInfo>& monomers)
{
    std::vector<std::string> symbols;
    for (const auto& info : monomers) {
        symbols.push_back(info.symbol.value_or("") + "/" +
                          info.polymer_type.value_or(""));
    }
    return symbols;
}

/**
 * Check that the filter's effective settings are the most permissive ones for
 * the monomer type and natural analog, i.e. only the polymer types are applied
 */
static void
check_only_polymer_types_applied(const MonomerDatabaseFilter& filter)
{
    BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::ANY);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "");
    BOOST_TEST(!filter.isMonomerTypeFiltered());
    BOOST_TEST(!filter.isNaturalAnalogFiltered());
    BOOST_TEST(filter.getNumFiltersApplied() ==
               (filter.isPolymerTypeFiltered() ? 1 : 0));
}

/**
 * Check that the filter applies {Peptide} + Backbone + C
 */
static void check_peptide_backbone_cys(const MonomerDatabaseFilter& filter)
{
    BOOST_TEST(filter.getPolymerTypes() == PEPTIDE);
    BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::BACKBONE);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::SPECIFIC);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "C");
    BOOST_TEST(filter.getNumFiltersApplied() == 3);
    BOOST_TEST(filter.getDescription() == "Peptide, backbone, Cys analogs");
}

BOOST_AUTO_TEST_CASE(test_polymer_type_flags)
{
    BOOST_TEST(has_flag(PolymerType::ALL, PEPTIDE));
    BOOST_TEST(has_flag(PEPTIDE | CHEM, CHEM));
    BOOST_TEST(!has_flag(PEPTIDE | CHEM, NUCLEIC_ACID));
    BOOST_TEST(!has_flag(PolymerType::ALL, PolymerType::NONE));
    BOOST_TEST(!has_flag(PEPTIDE, PEPTIDE | CHEM));
    BOOST_TEST(count_flags(PolymerType::NONE) == 0);
    BOOST_TEST(count_flags(NUCLEIC_ACID) == 1);
    BOOST_TEST(count_flags(PEPTIDE | NUCLEIC_ACID) == 2);
    BOOST_TEST(count_flags(PolymerType::ALL) == 3);
    BOOST_TEST((PEPTIDE | NUCLEIC_ACID | CHEM) == PolymerType::ALL);
    BOOST_TEST(((PEPTIDE | CHEM) & (NUCLEIC_ACID | CHEM)) == CHEM);
}

BOOST_AUTO_TEST_CASE(test_get_monomer_db_polymer_type)
{
    BOOST_TEST(get_monomer_db_polymer_type(make_monomer("A", "PEPTIDE")) ==
               PEPTIDE);
    BOOST_TEST(get_monomer_db_polymer_type(make_monomer("A", "RNA")) ==
               NUCLEIC_ACID);
    BOOST_TEST(get_monomer_db_polymer_type(make_monomer("A", "CHEM")) == CHEM);
    BOOST_TEST(get_monomer_db_polymer_type(make_monomer("A", "BOGUS")) ==
               PolymerType::NONE);
    BOOST_TEST(get_monomer_db_polymer_type(make_monomer("A", std::nullopt)) ==
               PolymerType::NONE);
    BOOST_TEST(get_polymer_type_display_text(PEPTIDE) == "Peptide");
    BOOST_TEST(get_polymer_type_display_text(NUCLEIC_ACID) == "Nucleic acid");
    BOOST_TEST(get_polymer_type_display_text(CHEM) == "Structure");
    BOOST_TEST(get_polymer_type_display_text(PolymerType::NONE) == "");
}

BOOST_AUTO_TEST_CASE(test_get_monomer_database_table_monomers)
{
    std::vector<MonomerInfo> custom{
        make_monomer("dC", "PEPTIDE", "C"),
        make_monomer("5mC", "RNA", "C"),
        make_monomer("Eth", "CHEM"),
        // replaces the core "ac"
        make_monomer("ac", "PEPTIDE", "C", std::nullopt, "custom ac"),
        // doesn't replace the core "am" because the polymer type differs
        make_monomer("am", "RNA", "A"),
    };
    std::vector<MonomerInfo> core{
        make_monomer("A", "PEPTIDE", "A"),
        make_monomer("ac", "PEPTIDE", "X", std::nullopt, "core ac"),
        make_monomer("am", "PEPTIDE", "X"),
        make_monomer("noAnalog", "PEPTIDE", ""),
        make_monomer("missingAnalog", "PEPTIDE"),
        // N is asparagine for peptides, so this has a natural analog
        make_monomer("dN", "PEPTIDE", "N"),
        make_monomer("N", "RNA", "N"),
        // X is a natural analog for nucleic acids, so this isn't excluded
        // as "None"
        make_monomer("rX", "RNA", "X"),
        make_monomer("R", "RNA", "R"),
        make_monomer("coreChem", "CHEM", "R"),
    };

    auto table_monomers = get_monomer_database_table_monomers(custom, core);
    std::vector<std::string> expected{
        "dC/PEPTIDE",
        "5mC/RNA",
        "Eth/CHEM",
        "ac/PEPTIDE",
        "am/RNA",
        "am/PEPTIDE",
        "noAnalog/PEPTIDE",
        "missingAnalog/PEPTIDE",
        "N/RNA",
        "coreChem/CHEM",
    };
    BOOST_TEST(get_symbols(table_monomers) == expected,
               boost::test_tools::per_element());
    // the custom "ac" is kept, not the core one
    BOOST_TEST(table_monomers[3].name.value_or("") == "custom ac");

    // with no custom monomers, only the core monomers with no natural analog
    // are included
    table_monomers = get_monomer_database_table_monomers({}, core);
    expected = {"ac/PEPTIDE",
                "am/PEPTIDE",
                "noAnalog/PEPTIDE",
                "missingAnalog/PEPTIDE",
                "N/RNA",
                "coreChem/CHEM"};
    BOOST_TEST(get_symbols(table_monomers) == expected,
               boost::test_tools::per_element());
}

BOOST_DATA_TEST_CASE(
    test_polymer_type_matching,
    boost::unit_test::data::make(ALL_POLYMER_TYPE_COMBINATIONS), polymer_types)
{
    MonomerDatabaseFilter filter(polymer_types);
    BOOST_TEST(filter.getPolymerTypes() == polymer_types);
    BOOST_TEST(filter.matches(make_monomer("A", "PEPTIDE", "A")) ==
               has_flag(polymer_types, PEPTIDE));
    BOOST_TEST(filter.matches(make_monomer("A", "RNA", "A")) ==
               has_flag(polymer_types, NUCLEIC_ACID));
    BOOST_TEST(filter.matches(make_monomer("Eth", "CHEM")) ==
               has_flag(polymer_types, CHEM));
    // rows with a missing or unparseable polymer type never match
    BOOST_TEST(!filter.matches(make_monomer("A", "BOGUS", "A")));
    BOOST_TEST(!filter.matches(make_monomer("A", std::nullopt, "A")));
    BOOST_TEST(filter.isPolymerTypeFiltered() ==
               (polymer_types != PolymerType::ALL));
}

BOOST_DATA_TEST_CASE(
    test_applicability,
    boost::unit_test::data::make(ALL_POLYMER_TYPE_COMBINATIONS), polymer_types)
{
    bool peptide = has_flag(polymer_types, PEPTIDE);
    bool nucleic_acid = has_flag(polymer_types, NUCLEIC_ACID);
    bool monomer_type_applicable = peptide && !nucleic_acid;
    bool analog_applicable = peptide != nucleic_acid;

    // start from {Peptide} + Backbone + C, then switch to the polymer types
    MonomerDatabaseFilter filter(PEPTIDE);
    filter.setMonomerType(MonomerTypeFilter::BACKBONE);
    filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    filter.setPolymerTypes(polymer_types);

    BOOST_TEST(filter.isMonomerTypeApplicable() == monomer_type_applicable);
    BOOST_TEST(filter.isNaturalAnalogApplicable() == analog_applicable);
    BOOST_TEST(filter.getMonomerType() == (monomer_type_applicable
                                               ? MonomerTypeFilter::BACKBONE
                                               : MonomerTypeFilter::ANY));
    if (analog_applicable && peptide) {
        BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::SPECIFIC);
        BOOST_TEST(filter.getSpecificNaturalAnalog() == "C");
    } else {
        // either inapplicable, or applicable for nucleic acids, in which case
        // the peptide analog is forgotten
        BOOST_TEST(filter.getNaturalAnalog() ==
                   NaturalAnalogFilter::ANY_OR_NONE);
        BOOST_TEST(filter.getSpecificNaturalAnalog() == "");
    }
    BOOST_TEST(filter.isMonomerTypeFiltered() == monomer_type_applicable);
    BOOST_TEST(filter.isNaturalAnalogFiltered() ==
               (analog_applicable && peptide));
}

BOOST_AUTO_TEST_CASE(test_natural_analog_choices)
{
    MonomerDatabaseFilter filter(PEPTIDE | CHEM);
    std::vector<std::string> expected{"A", "R", "N", "D", "C", "Q", "E",
                                      "G", "H", "I", "L", "K", "M", "F",
                                      "P", "S", "T", "W", "Y", "V"};
    BOOST_TEST(filter.getNaturalAnalogChoices() == expected,
               boost::test_tools::per_element());

    filter.setPolymerTypes(NUCLEIC_ACID);
    expected = {"R", "P", "A", "C", "G", "U", "T"};
    BOOST_TEST(filter.getNaturalAnalogChoices() == expected,
               boost::test_tools::per_element());

    for (auto polymer_types :
         {PolymerType::NONE, CHEM, PEPTIDE | NUCLEIC_ACID}) {
        filter.setPolymerTypes(polymer_types);
        BOOST_TEST(filter.getNaturalAnalogChoices().empty());
    }
}

/**
 * {Peptide} + Backbone + C, then {Peptide, Nucleic acid}, then {Peptide}
 * restores Backbone and C. The inapplicable middle state matches, counts and
 * describes as Any / Any or None.
 */
BOOST_AUTO_TEST_CASE(test_remembered_through_peptide_and_nucleic_acid)
{
    MonomerDatabaseFilter filter(PEPTIDE);
    filter.setMonomerType(MonomerTypeFilter::BACKBONE);
    filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    check_peptide_backbone_cys(filter);
    auto terminal_ala = make_monomer("am", "PEPTIDE", "A", "Terminal");
    BOOST_TEST(!filter.matches(terminal_ala));

    filter.setPolymerTypes(PEPTIDE | NUCLEIC_ACID);
    check_only_polymer_types_applied(filter);
    BOOST_TEST(filter.getNumFiltersApplied() == 1);
    BOOST_TEST(filter.getDescription() == "Peptide or nucleic acid monomers");
    BOOST_TEST(filter.matches(terminal_ala));
    BOOST_TEST(filter.matches(make_monomer("A", "RNA", "A")));

    filter.setPolymerTypes(PEPTIDE);
    check_peptide_backbone_cys(filter);
    BOOST_TEST(!filter.matches(terminal_ala));
}

/**
 * The same as above, through {Peptide, Structure} -> {Structure} ->
 * {Peptide, Structure}
 */
BOOST_AUTO_TEST_CASE(test_remembered_through_structure)
{
    MonomerDatabaseFilter filter(PEPTIDE | CHEM);
    filter.setMonomerType(MonomerTypeFilter::BACKBONE);
    filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    BOOST_TEST(filter.getDescription() == "Peptide, backbone, Cys analogs");
    auto chem = make_monomer("Eth", "CHEM");
    BOOST_TEST(!filter.matches(chem));

    filter.setPolymerTypes(CHEM);
    check_only_polymer_types_applied(filter);
    BOOST_TEST(filter.getDescription() == "Structure monomers");
    BOOST_TEST(filter.matches(chem));

    filter.setPolymerTypes(PEPTIDE | CHEM);
    BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::BACKBONE);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::SPECIFIC);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "C");
    BOOST_TEST(filter.getDescription() == "Peptide, backbone, Cys analogs");
    BOOST_TEST(!filter.matches(chem));
}

/**
 * A peptide analog (C, None or Any) is forgotten when the analog applies again
 * for Nucleic acid, through either toggle order, and stays forgotten on
 * returning to {Peptide}
 */
BOOST_DATA_TEST_CASE(test_peptide_analog_forgotten,
                     boost::unit_test::data::make(
                         {NaturalAnalogFilter::SPECIFIC,
                          NaturalAnalogFilter::NONE, NaturalAnalogFilter::ANY}),
                     analog)
{
    // {Peptide} -> {Peptide, Nucleic acid} -> {Nucleic acid}
    MonomerDatabaseFilter filter(PEPTIDE);
    filter.setNaturalAnalog(analog, "C");
    filter.setPolymerTypes(PEPTIDE | NUCLEIC_ACID);
    filter.setPolymerTypes(NUCLEIC_ACID);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    filter.setPolymerTypes(PEPTIDE);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "");
    BOOST_TEST((filter == MonomerDatabaseFilter(PEPTIDE)));

    // {Peptide, Structure} -> {Structure} -> {Nucleic acid, Structure}
    filter = MonomerDatabaseFilter(PEPTIDE | CHEM);
    filter.setNaturalAnalog(analog, "C");
    filter.setPolymerTypes(CHEM);
    filter.setPolymerTypes(NUCLEIC_ACID | CHEM);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    filter.setPolymerTypes(PEPTIDE);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "");
}

/**
 * The reverse of the above: a nucleic-acid analog is forgotten when the analog
 * applies again for Peptide, but restored when it applies again for Nucleic
 * acid
 */
BOOST_DATA_TEST_CASE(test_nucleic_acid_analog_forgotten,
                     boost::unit_test::data::make(
                         {NaturalAnalogFilter::SPECIFIC,
                          NaturalAnalogFilter::NONE, NaturalAnalogFilter::ANY}),
                     analog)
{
    MonomerDatabaseFilter filter(NUCLEIC_ACID);
    filter.setNaturalAnalog(analog, "A");
    filter.setPolymerTypes(PEPTIDE | NUCLEIC_ACID);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    // returning to the same polymer type restores the analog
    filter.setPolymerTypes(NUCLEIC_ACID);
    BOOST_TEST(filter.getNaturalAnalog() == analog);
    BOOST_TEST(filter.getSpecificNaturalAnalog() ==
               (analog == NaturalAnalogFilter::SPECIFIC ? "A" : ""));

    // {Nucleic acid} -> {Peptide, Nucleic acid} -> {Peptide}
    filter.setPolymerTypes(PEPTIDE | NUCLEIC_ACID);
    filter.setPolymerTypes(PEPTIDE);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    filter.setPolymerTypes(NUCLEIC_ACID);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);

    // {Nucleic acid, Structure} -> {Structure} -> {Peptide, Structure}
    filter = MonomerDatabaseFilter(NUCLEIC_ACID | CHEM);
    filter.setNaturalAnalog(analog, "A");
    filter.setPolymerTypes(CHEM);
    filter.setPolymerTypes(PEPTIDE | CHEM);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
    filter.setPolymerTypes(NUCLEIC_ACID | CHEM);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::ANY_OR_NONE);
}

/**
 * Monomer type survives {Peptide, Nucleic acid} and {Structure}, but is
 * forgotten on reaching {Nucleic acid} or {Nucleic acid, Structure}, and stays
 * Any on returning to {Peptide}
 */
BOOST_AUTO_TEST_CASE(test_monomer_type_remembered_and_forgotten)
{
    for (auto intermediate : {PEPTIDE | NUCLEIC_ACID, CHEM, PolymerType::ALL}) {
        MonomerDatabaseFilter filter(PEPTIDE);
        filter.setMonomerType(MonomerTypeFilter::TERMINAL);
        filter.setPolymerTypes(intermediate);
        BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::ANY);
        filter.setPolymerTypes(PEPTIDE);
        BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::TERMINAL);
    }

    for (auto path : std::vector<std::vector<PolymerType>>{
             {NUCLEIC_ACID},
             {PEPTIDE | NUCLEIC_ACID, NUCLEIC_ACID},
             {CHEM, NUCLEIC_ACID | CHEM},
             {PEPTIDE | NUCLEIC_ACID | CHEM, NUCLEIC_ACID | CHEM}}) {
        MonomerDatabaseFilter filter(PEPTIDE | CHEM);
        filter.setMonomerType(MonomerTypeFilter::TERMINAL);
        for (auto polymer_types : path) {
            filter.setPolymerTypes(polymer_types);
            BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::ANY);
        }
        filter.setPolymerTypes(PEPTIDE);
        BOOST_TEST(filter.getMonomerType() == MonomerTypeFilter::ANY);
        BOOST_TEST((filter == MonomerDatabaseFilter(PEPTIDE)));
    }
}

/**
 * setNaturalAnalog records the polymer type that the analog was applicable
 * for, ignoring Structure
 */
BOOST_AUTO_TEST_CASE(test_set_natural_analog_records_polymer_type)
{
    // recorded for Peptide even though Structure was also checked
    MonomerDatabaseFilter filter(PEPTIDE | CHEM);
    filter.setNaturalAnalog(NaturalAnalogFilter::ANY_OR_NONE);
    BOOST_TEST(!(filter == MonomerDatabaseFilter(PEPTIDE | CHEM)));
    filter.setPolymerTypes(PEPTIDE);
    BOOST_TEST(!(filter == MonomerDatabaseFilter(PEPTIDE)));
    // the recorded type is forgotten along with the value
    filter.setPolymerTypes(NUCLEIC_ACID);
    BOOST_TEST((filter == MonomerDatabaseFilter(NUCLEIC_ACID)));

    // recorded for Nucleic acid, so it survives a round trip through
    // {Structure}
    filter = MonomerDatabaseFilter(NUCLEIC_ACID | CHEM);
    filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "G");
    filter.setPolymerTypes(CHEM);
    filter.setPolymerTypes(NUCLEIC_ACID);
    BOOST_TEST(filter.getNaturalAnalog() == NaturalAnalogFilter::SPECIFIC);
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "G");

    // the specific analog is only stored for SPECIFIC
    filter.setNaturalAnalog(NaturalAnalogFilter::NONE, "G");
    BOOST_TEST(filter.getSpecificNaturalAnalog() == "");
    MonomerDatabaseFilter other(NUCLEIC_ACID);
    other.setNaturalAnalog(NaturalAnalogFilter::NONE);
    BOOST_TEST((filter == other));
}

/**
 * == sees remembered values, and the effective getters don't
 */
BOOST_AUTO_TEST_CASE(test_equality)
{
    BOOST_TEST(
        (MonomerDatabaseFilter() == MonomerDatabaseFilter(PolymerType::ALL)));
    BOOST_TEST(!(MonomerDatabaseFilter() == MonomerDatabaseFilter(PEPTIDE)));

    MonomerDatabaseFilter remembering(PEPTIDE);
    remembering.setMonomerType(MonomerTypeFilter::BACKBONE);
    remembering.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    remembering.setPolymerTypes(PEPTIDE | NUCLEIC_ACID);
    MonomerDatabaseFilter plain(PEPTIDE | NUCLEIC_ACID);

    BOOST_TEST(remembering.getPolymerTypes() == plain.getPolymerTypes());
    BOOST_TEST(remembering.getMonomerType() == plain.getMonomerType());
    BOOST_TEST(remembering.getNaturalAnalog() == plain.getNaturalAnalog());
    BOOST_TEST(remembering.getSpecificNaturalAnalog() ==
               plain.getSpecificNaturalAnalog());
    BOOST_TEST(remembering.getDescription() == plain.getDescription());
    BOOST_TEST(!(remembering == plain));

    MonomerDatabaseFilter expected(PEPTIDE);
    expected.setMonomerType(MonomerTypeFilter::BACKBONE);
    expected.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    remembering.setPolymerTypes(PEPTIDE);
    BOOST_TEST((remembering == expected));
}

BOOST_AUTO_TEST_CASE(test_chem_rows)
{
    auto chem = make_monomer("Eth", "CHEM", "X", "Undefined");
    auto chem_no_analog = make_monomer("Eth2", "CHEM");

    MonomerDatabaseFilter filter(PEPTIDE | CHEM);
    BOOST_TEST(filter.matches(chem));
    BOOST_TEST(filter.matches(chem_no_analog));
    for (auto monomer_type :
         {MonomerTypeFilter::BACKBONE, MonomerTypeFilter::TERMINAL}) {
        filter.setMonomerType(monomer_type);
        BOOST_TEST(!filter.matches(chem));
        BOOST_TEST(!filter.matches(chem_no_analog));
    }
    filter.setMonomerType(MonomerTypeFilter::ANY);

    for (auto polymer_types : {PEPTIDE | CHEM, NUCLEIC_ACID | CHEM}) {
        filter.setPolymerTypes(polymer_types);
        for (auto [analog, expected] :
             std::vector<std::pair<NaturalAnalogFilter, bool>>{
                 {NaturalAnalogFilter::ANY_OR_NONE, true},
                 {NaturalAnalogFilter::NONE, true},
                 {NaturalAnalogFilter::ANY, false},
                 {NaturalAnalogFilter::SPECIFIC, false},
             }) {
            filter.setNaturalAnalog(analog, "X");
            BOOST_TEST(filter.matches(chem) == expected);
            BOOST_TEST(filter.matches(chem_no_analog) == expected);
        }
    }
    BOOST_TEST(has_no_natural_analog(chem));
    BOOST_TEST(get_natural_analog_display_text(chem) == "");
}

BOOST_AUTO_TEST_CASE(test_none_and_any_analogs)
{
    struct AnalogCase {
        MonomerInfo info;
        PolymerType polymer_type;
        bool has_no_analog;
    };
    std::vector<AnalogCase> cases{
        {make_monomer("pX", "PEPTIDE", "X"), PEPTIDE, true},
        {make_monomer("pEmpty", "PEPTIDE", ""), PEPTIDE, true},
        {make_monomer("pMissing", "PEPTIDE"), PEPTIDE, true},
        // N is asparagine for peptides
        {make_monomer("pN", "PEPTIDE", "N"), PEPTIDE, false},
        {make_monomer("pC", "PEPTIDE", "C"), PEPTIDE, false},
        {make_monomer("nN", "RNA", "N"), NUCLEIC_ACID, true},
        {make_monomer("nEmpty", "RNA", ""), NUCLEIC_ACID, true},
        {make_monomer("nMissing", "RNA"), NUCLEIC_ACID, true},
        {make_monomer("nX", "RNA", "X"), NUCLEIC_ACID, false},
        {make_monomer("nA", "RNA", "A"), NUCLEIC_ACID, false},
    };
    for (const auto& [info, polymer_type, has_no_analog] : cases) {
        BOOST_TEST_CONTEXT(*info.symbol)
        {
            BOOST_TEST(has_no_natural_analog(info) == has_no_analog);
            BOOST_TEST(get_natural_analog_display_text(info) ==
                       (has_no_analog ? "" : *info.natural_analog));

            MonomerDatabaseFilter filter(polymer_type);
            BOOST_TEST(filter.matches(info));
            filter.setNaturalAnalog(NaturalAnalogFilter::NONE);
            BOOST_TEST(filter.matches(info) == has_no_analog);
            filter.setNaturalAnalog(NaturalAnalogFilter::ANY);
            BOOST_TEST(filter.matches(info) == !has_no_analog);
            filter.setNaturalAnalog(NaturalAnalogFilter::ANY_OR_NONE);
            BOOST_TEST(filter.matches(info));
            filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC,
                                    info.natural_analog.value_or(""));
            BOOST_TEST(filter.matches(info) == !has_no_analog);
        }
    }

    // a specific analog doesn't match other analogs
    MonomerDatabaseFilter filter(PEPTIDE);
    filter.setNaturalAnalog(NaturalAnalogFilter::SPECIFIC, "C");
    BOOST_TEST(filter.matches(make_monomer("dC", "PEPTIDE", "C")));
    BOOST_TEST(!filter.matches(make_monomer("dA", "PEPTIDE", "A")));
    BOOST_TEST(!filter.matches(make_monomer("c", "PEPTIDE", "c")));
}

BOOST_AUTO_TEST_CASE(test_monomer_type_case_insensitive)
{
    MonomerDatabaseFilter filter(PEPTIDE);
    filter.setMonomerType(MonomerTypeFilter::BACKBONE);
    for (auto monomer_type : {"Backbone", "backbone", "BACKBONE"}) {
        BOOST_TEST(
            filter.matches(make_monomer("A", "PEPTIDE", "A", monomer_type)));
    }
    BOOST_TEST(!filter.matches(make_monomer("A", "PEPTIDE", "A", "Terminal")));
    BOOST_TEST(!filter.matches(make_monomer("A", "PEPTIDE", "A", "Branch")));
    BOOST_TEST(!filter.matches(make_monomer("A", "PEPTIDE", "A")));

    filter.setMonomerType(MonomerTypeFilter::TERMINAL);
    for (auto monomer_type : {"Terminal", "terminal", "TERMINAL"}) {
        BOOST_TEST(
            filter.matches(make_monomer("A", "PEPTIDE", "A", monomer_type)));
    }
    BOOST_TEST(!filter.matches(make_monomer("A", "PEPTIDE", "A", "Backbone")));
}

/**
 * Monomers used for the sort and search tests. The comments give the values
 * that are sorted on.
 */
static std::vector<MonomerInfo> get_sort_monomers()
{
    return {
        // 0: symbol b, name Zeta, Peptide, analog C
        make_monomer("b", "PEPTIDE", "C", std::nullopt, "Zeta"),
        // 1: symbol A, name alpha, Nucleic acid, analog A
        make_monomer("A", "RNA", "A", std::nullopt, "alpha"),
        // 2: symbol c, no name, Structure, no analog
        make_monomer("c", "CHEM", std::nullopt, std::nullopt, std::nullopt),
        // 3: symbol a, name Beta, Peptide, no analog (X)
        make_monomer("a", "PEPTIDE", "X", std::nullopt, "Beta"),
        // 4: symbol A, name alpha, Peptide, analog A
        make_monomer("A", "PEPTIDE", "A", std::nullopt, "alpha"),
        // 5: no symbol, name Gamma, Peptide, analog G
        make_monomer(std::nullopt, "PEPTIDE", "G", std::nullopt, "Gamma"),
    };
}

struct SortCase {
    MonomerDbColumn column;
    bool ascending;
    std::vector<std::size_t> expected;
};

std::ostream& operator<<(std::ostream& os, const SortCase& sort_case)
{
    os << "column " << static_cast<int>(sort_case.column)
       << (sort_case.ascending ? " ascending" : " descending");
    return os;
}

// Ties are broken by symbol and then by polymer type ("Nucleic acid" before
// "Peptide"). Monomers 3 and 4 tie on everything (case-insensitively), so they
// keep their original order. Missing values sort last when ascending and first
// when descending, and tie-breaks are always ascending.
static const std::vector<SortCase> SORT_CASES{
    {MonomerDbColumn::SYMBOL, true, {1, 3, 4, 0, 2, 5}},
    {MonomerDbColumn::SYMBOL, false, {5, 2, 0, 1, 3, 4}},
    {MonomerDbColumn::NAME, true, {1, 4, 3, 5, 0, 2}},
    {MonomerDbColumn::NAME, false, {2, 0, 5, 3, 1, 4}},
    {MonomerDbColumn::POLYMER_TYPE, true, {1, 3, 4, 0, 5, 2}},
    {MonomerDbColumn::POLYMER_TYPE, false, {2, 3, 4, 0, 5, 1}},
    // None analogs (X and missing) are blank, so they sort as missing values
    {MonomerDbColumn::NATURAL_ANALOG, true, {1, 4, 0, 5, 3, 2}},
    {MonomerDbColumn::NATURAL_ANALOG, false, {3, 2, 5, 0, 1, 4}},
};

BOOST_DATA_TEST_CASE(test_sort, boost::unit_test::data::make(SORT_CASES),
                     sort_case)
{
    auto monomers = get_sort_monomers();
    auto rows = filter_and_sort_monomers(
        monomers, MonomerDatabaseFilter(),
        MonomerDbSort{sort_case.column, sort_case.ascending});
    BOOST_TEST(rows == sort_case.expected, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(test_filter_and_sort)
{
    auto monomers = get_sort_monomers();
    auto rows = filter_and_sort_monomers(
        monomers, MonomerDatabaseFilter(PEPTIDE), MonomerDbSort{});
    std::vector<std::size_t> expected{3, 4, 0, 5};
    BOOST_TEST(rows == expected, boost::test_tools::per_element());

    MonomerDatabaseFilter filter(PEPTIDE | CHEM);
    filter.setNaturalAnalog(NaturalAnalogFilter::NONE);
    rows = filter_and_sort_monomers(
        monomers, filter, MonomerDbSort{MonomerDbColumn::SYMBOL, false});
    expected = {2, 3};
    BOOST_TEST(rows == expected, boost::test_tools::per_element());

    BOOST_TEST(filter_and_sort_monomers({}, filter, MonomerDbSort{}).empty());
}

BOOST_AUTO_TEST_CASE(test_search)
{
    auto monomers = get_sort_monomers();
    std::vector<std::size_t> rows{5, 3, 0, 1, 2, 4};

    // empty text returns the input unchanged
    BOOST_TEST(filter_by_search_text(monomers, rows, "") == rows,
               boost::test_tools::per_element());
    // name only, case-insensitively, preserving order
    std::vector<std::size_t> expected{3, 0};
    BOOST_TEST(filter_by_search_text(monomers, rows, "ETA") == expected,
               boost::test_tools::per_element());
    // symbol only (monomer 2 has no name)
    expected = {2};
    BOOST_TEST(filter_by_search_text(monomers, rows, "C") == expected,
               boost::test_tools::per_element());
    // symbol or name
    expected = {3, 0};
    BOOST_TEST(filter_by_search_text(monomers, rows, "b") == expected,
               boost::test_tools::per_element());
    // only the given rows are searched
    expected = {4};
    BOOST_TEST(filter_by_search_text(monomers, {2, 4}, "alpha") == expected,
               boost::test_tools::per_element());
    BOOST_TEST(filter_by_search_text(monomers, rows, "nothing").empty());

    // substring search, not prefix search
    std::vector<MonomerInfo> hydro_monomers{
        make_monomer("Hyl", "PEPTIDE", "K", std::nullopt, "Hydroxylysine"),
        make_monomer("Ile(4OH)", "PEPTIDE", "I", std::nullopt,
                     "Hydroxyisoleucine"),
        make_monomer("Me_dPhe", "PEPTIDE", "F", std::nullopt,
                     "Me_Phe(a,b-dehydro)"),
        make_monomer("A", "PEPTIDE", "A", std::nullopt, "Alanine"),
    };
    expected = {0, 1, 2};
    BOOST_TEST(filter_by_search_text(hydro_monomers, {0, 1, 2, 3}, "hydro") ==
                   expected,
               boost::test_tools::per_element());
}

struct DescriptionCase {
    PolymerType polymer_types;
    MonomerTypeFilter monomer_type;
    NaturalAnalogFilter analog;
    std::string specific;
    std::string description;
};

std::ostream& operator<<(std::ostream& os, const DescriptionCase& desc_case)
{
    os << desc_case.description;
    return os;
}

static const std::vector<DescriptionCase> DESCRIPTION_CASES{
    {PEPTIDE, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY_OR_NONE, "",
     "Peptide monomers"},
    {NUCLEIC_ACID, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY_OR_NONE, "",
     "Nucleic acid monomers"},
    {CHEM, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY_OR_NONE, "",
     "Structure monomers"},
    {PEPTIDE | CHEM, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY_OR_NONE,
     "", "Peptide or structure monomers"},
    {NUCLEIC_ACID | CHEM, MonomerTypeFilter::ANY,
     NaturalAnalogFilter::ANY_OR_NONE, "",
     "Nucleic acid or structure monomers"},
    {PEPTIDE | NUCLEIC_ACID, MonomerTypeFilter::ANY,
     NaturalAnalogFilter::ANY_OR_NONE, "", "Peptide or nucleic acid monomers"},
    {PolymerType::ALL, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY_OR_NONE,
     "", "All monomers"},
    {PEPTIDE, MonomerTypeFilter::TERMINAL, NaturalAnalogFilter::ANY_OR_NONE, "",
     "Peptide, terminal monomers"},
    {PEPTIDE, MonomerTypeFilter::ANY, NaturalAnalogFilter::SPECIFIC, "C",
     "Peptide, Cys analogs"},
    {PEPTIDE, MonomerTypeFilter::BACKBONE, NaturalAnalogFilter::SPECIFIC, "C",
     "Peptide, backbone, Cys analogs"},
    {NUCLEIC_ACID, MonomerTypeFilter::ANY, NaturalAnalogFilter::SPECIFIC, "A",
     "Nucleic acid, A analogs"},
    {PEPTIDE | CHEM, MonomerTypeFilter::TERMINAL,
     NaturalAnalogFilter::ANY_OR_NONE, "", "Peptide, terminal monomers"},
    {PEPTIDE | CHEM, MonomerTypeFilter::ANY, NaturalAnalogFilter::SPECIFIC, "C",
     "Peptide, Cys analogs"},
    {PEPTIDE, MonomerTypeFilter::ANY, NaturalAnalogFilter::NONE, "",
     "Peptide monomers with no natural analog"},
    {PEPTIDE | CHEM, MonomerTypeFilter::ANY, NaturalAnalogFilter::NONE, "",
     "Peptide or structure monomers with no natural analog"},
    {PEPTIDE, MonomerTypeFilter::TERMINAL, NaturalAnalogFilter::NONE, "",
     "Peptide, terminal monomers with no natural analog"},
    {PEPTIDE, MonomerTypeFilter::ANY, NaturalAnalogFilter::ANY, "",
     "Peptide monomers with any natural analog"},
};

BOOST_DATA_TEST_CASE(test_count_and_description,
                     boost::unit_test::data::make(DESCRIPTION_CASES), desc_case)
{
    MonomerDatabaseFilter filter(desc_case.polymer_types);
    filter.setMonomerType(desc_case.monomer_type);
    if (filter.isNaturalAnalogApplicable()) {
        filter.setNaturalAnalog(desc_case.analog, desc_case.specific);
    }
    BOOST_TEST(filter.getDescription() == desc_case.description);

    bool polymer_type_filtered = desc_case.polymer_types != PolymerType::ALL;
    bool monomer_type_filtered =
        desc_case.monomer_type != MonomerTypeFilter::ANY;
    bool analog_filtered = desc_case.analog != NaturalAnalogFilter::ANY_OR_NONE;
    int count = polymer_type_filtered + monomer_type_filtered + analog_filtered;
    BOOST_TEST(filter.isPolymerTypeFiltered() == polymer_type_filtered);
    BOOST_TEST(filter.isMonomerTypeFiltered() == monomer_type_filtered);
    BOOST_TEST(filter.isNaturalAnalogFiltered() == analog_filtered);
    BOOST_TEST(filter.getNumFiltersApplied() == count);
    BOOST_TEST(filter.isFiltering() == (count > 0));
}

BOOST_AUTO_TEST_CASE(test_default_filters)
{
    // the most permissive filter (Clear All)
    MonomerDatabaseFilter most_permissive;
    BOOST_TEST(most_permissive.getPolymerTypes() == PolymerType::ALL);
    BOOST_TEST(most_permissive.getNumFiltersApplied() == 0);
    BOOST_TEST(!most_permissive.isFiltering());
    BOOST_TEST(!most_permissive.isPolymerTypeFiltered());
    BOOST_TEST(most_permissive.getDescription() == "All monomers");

    // the page defaults
    for (auto [polymer_types, description] :
         std::vector<std::pair<PolymerType, std::string>>{
             {PEPTIDE | CHEM, "Peptide or structure monomers"},
             {NUCLEIC_ACID | CHEM, "Nucleic acid or structure monomers"}}) {
        MonomerDatabaseFilter page_default(polymer_types);
        BOOST_TEST(page_default.getNumFiltersApplied() == 1);
        BOOST_TEST(page_default.isFiltering());
        BOOST_TEST(page_default.isPolymerTypeFiltered());
        BOOST_TEST(!page_default.isMonomerTypeFiltered());
        BOOST_TEST(!page_default.isNaturalAnalogFiltered());
        BOOST_TEST(page_default.isNaturalAnalogApplicable());
        BOOST_TEST(page_default.isMonomerTypeApplicable() ==
                   (polymer_types == (PEPTIDE | CHEM)));
        BOOST_TEST(page_default.getDescription() == description);
    }
}

} // namespace sketcher
} // namespace schrodinger
