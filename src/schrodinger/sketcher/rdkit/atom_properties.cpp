#include "schrodinger/sketcher/rdkit/atom_properties.h"

#include <deque>
#include <unordered_set>
#include <initializer_list>

#include <boost/algorithm/string.hpp>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"

namespace schrodinger::sketcher
{

// In AtomType queries, atomic numbers above this value indicate that the atom
// is aromatic
const int ATOM_TYPE_AROMATIC_OFFSET = 1000;

const std::unordered_map<std::string, AtomQuery> TYPE_LABEL_TO_ATOM_QUERY = {
    {"A", AtomQuery::A},   {"AH", AtomQuery::AH}, {"Q", AtomQuery::Q},
    {"QH", AtomQuery::QH}, {"M", AtomQuery::M},   {"MH", AtomQuery::MH},
    {"X", AtomQuery::X},   {"XH", AtomQuery::XH},
};

// some properties are affected by more than one query type (e.g.
// AtomIsAliphatic and AtomIsAromatic both affect
// AtomQueryProperties::aromaticity). For those scenarios, we use one of
// these constants as the seen_descriptions key instead of the query
// description.
static const std::string ATOM_AROMATICITY = "ATOM_AROMATICITY";
static const std::string ATOM_RING_COUNT = "ATOM_RING_COUNT";
static const std::string QUERY_TYPE = "QUERY_TYPE";

using AtomQueryPropertyPtr = std::variant<
    unsigned int AtomQueryProperties::*,
    std::optional<unsigned int> AtomQueryProperties::*,
    std::optional<int> AtomQueryProperties::*,
    std::string AtomQueryProperties::*, Element AtomQueryProperties::*,
    std::vector<Element> AtomQueryProperties::*,
    std::optional<EnhancedStereo> AtomQueryProperties::*,
    AtomQuery AtomQueryProperties::*, QueryAromaticity AtomQueryProperties::*,
    QueryCount AtomQueryProperties::*, QueryType AtomQueryProperties::*>;
using AtomQueryPropertyList = std::vector<AtomQueryPropertyPtr>;

/**
 * Atom query properties that directly relate to the query type, i.e., the
 * properties from query type combo box itself, as well as the line edits and
 * spin boxs that appear/disappear depending on the query type selected in the
 * combo box.
 */
const AtomQueryPropertyList QUERY_TYPE_PROPERTIES = {
    &AtomQueryProperties::query_type,   &AtomQueryProperties::element,
    &AtomQueryProperties::allowed_list, &AtomQueryProperties::wildcard,
    &AtomQueryProperties::r_group,      &AtomQueryProperties::smarts_query,
};

/**
 * Atom query properties that are set on the General tab of the Edit Atom
 * Properties dialog, but that don't directly relate to the query type
 */
const AtomQueryPropertyList GENERAL_QUERY_PROPERTIES = {
    &AtomQueryProperties::isotope,
    &AtomQueryProperties::charge,
    &AtomQueryProperties::unpaired_electrons,
    &AtomQueryProperties::enhanced_stereo,
};

/**
 * Atom query properties that are set on the Advanced tab of the Edit Atom
 * Properties dialog
 */
const AtomQueryPropertyList ADVANCED_QUERY_PROPERTIES = {
    &AtomQueryProperties::total_h_type,
    &AtomQueryProperties::total_h_exact_val,
    &AtomQueryProperties::num_connections,
    &AtomQueryProperties::aromaticity,
    &AtomQueryProperties::ring_count_type,
    &AtomQueryProperties::ring_count_exact_val,
    &AtomQueryProperties::ring_bond_count_type,
    &AtomQueryProperties::ring_bond_count_exact_val,
    &AtomQueryProperties::smallest_ring_size,
};

/**
 * Concatenate multiple AtomQueryPropertyLists into a single vector
 */
static AtomQueryPropertyList
concat_query_props(std::initializer_list<const AtomQueryPropertyList> vectors)
{
    AtomQueryPropertyList concatenated;
    size_t vec_size = std::accumulate(
        vectors.begin(), vectors.end(), 0,
        [](auto current_sum, auto vec) { return current_sum + vec.size(); });
    concatenated.reserve(vec_size);
    for (auto& cur_vector : vectors) {
        concatenated.insert(concatenated.end(), cur_vector.begin(),
                            cur_vector.end());
    }
    return concatenated;
}

const AtomQueryPropertyList ALL_QUERY_PROPERTIES =
    concat_query_props({QUERY_TYPE_PROPERTIES, GENERAL_QUERY_PROPERTIES,
                        ADVANCED_QUERY_PROPERTIES});

/**
 * An error that indicates that the query being parsed cannot be represented in
 * the property fields of the Edit Atom Properties dialog.  In this scenario,
 * the query will instead be loaded into the dialog as a SMARTS string (with the
 * query type set to SMARTS).
 */
class UnrecognizedQueryError : public std::runtime_error
{
  public:
    using std::runtime_error::runtime_error;
};

EnhancedStereo::EnhancedStereo(rdkit_extensions::EnhancedStereo enh_stereo) :
    EnhancedStereo(enh_stereo.first, enh_stereo.second)
{
}

RDKit::StereoGroupType EnhancedStereo::type() const
{
    return first;
}

void EnhancedStereo::setType(RDKit::StereoGroupType group_type)
{
    first = group_type;
}

unsigned int EnhancedStereo::groupId() const
{
    return type() == RDKit::StereoGroupType::STEREO_ABSOLUTE ? 0 : second;
}

void EnhancedStereo::setGroupId(unsigned int group_id)
{
    second = group_id;
}

/**
 * Determine whether the two query properties contain the same value for the
 * specified member variable
 */
template <typename T>
bool compare_property(const T AtomQueryProperties::*member_variable,
                      const AtomQueryProperties* this_query,
                      const AtomQueryProperties* other_query)
{
    return this_query->*member_variable == other_query->*member_variable;
}

/**
 * Assign the value for the specified member variable from one set of query
 * properties to the other
 */
template <typename T>
void assign_property(T AtomQueryProperties::*member_variable,
                     AtomQueryProperties* to_query,
                     const AtomQueryProperties* from_query)
{
    to_query->*member_variable = from_query->*member_variable;
}

/**
 * Determime whether teh two query properties contain the same values for all
 * specified member variables
 */
static bool compare_query_properties(const AtomQueryProperties* this_query,
                                     const AtomQueryProperties* other_query,
                                     const AtomQueryPropertyList properties)
{
    for (auto member_variable : properties) {
        if (!std::visit(
                [&this_query, &other_query](auto&& arg) {
                    return compare_property(arg, this_query, other_query);
                },
                member_variable)) {
            return false;
        }
    }
    return true;
}

/**
 * Return a new AtomQueryProperties object that contains some of the properties
 * from the input query.  All other property values will be set to their default
 * value.
 *
 * @param query The input query to take property values from
 * @param properties The properties to take from the input query
 */
AtomQueryProperties
get_subset_of_properties(const AtomQueryProperties* query,
                         const AtomQueryPropertyList properties)
{
    AtomQueryProperties new_query;
    for (auto member_variable : properties) {
        std::visit([&query, &new_query](
                       auto&& arg) { assign_property(arg, &new_query, query); },
                   member_variable);
    }
    return new_query;
}

bool AbstractAtomProperties::operator==(
    const AbstractAtomProperties& other) const
{
    if (isQuery() != other.isQuery()) {
        return false;
    } else if (isQuery()) {
        auto* this_query = static_cast<const AtomQueryProperties*>(this);
        auto* other_query = static_cast<const AtomQueryProperties*>(&other);
        return compare_query_properties(this_query, other_query,
                                        ALL_QUERY_PROPERTIES);
    } else {
        auto* this_atom = static_cast<const AtomProperties*>(this);
        auto* other_atom = static_cast<const AtomProperties*>(&other);
        return this_atom->element == other_atom->element &&
               this_atom->isotope == other_atom->isotope &&
               this_atom->enhanced_stereo == other_atom->enhanced_stereo &&
               this_atom->charge == other_atom->charge &&
               this_atom->unpaired_electrons == other_atom->unpaired_electrons;
    }
}

bool AbstractAtomProperties::operator!=(
    const AbstractAtomProperties& other) const
{
    return !(*this == other);
}

template <typename T> static std::string output_optional(std::optional<T> value)
{
    return value.has_value() ? std::to_string(*value) : "None";
}

std::ostream& operator<<(std::ostream& os, const AbstractAtomProperties& props)
{
    // TODO: enhanced stereo
    if (!props.isQuery()) {
        auto* atom_props = static_cast<const AtomProperties*>(&props);
        os << "AtomProperties:\n"
           << "\tElement: "
           << atomic_number_to_symbol(static_cast<int>(atom_props->element))
           << "\n"
           << "\tIsotope: " << output_optional(atom_props->isotope) << "\n"
           << "\tCharge: " << atom_props->charge << "\n"
           << "\tUnpaired electrons: " << atom_props->unpaired_electrons
           << "\n";
    } else {
        auto* query_props = static_cast<const AtomQueryProperties*>(&props);

        // convert the allowed list (or disallowed list) into a string
        std::string element_list_text;
        for (auto cur_elem : query_props->allowed_list) {
            if (cur_elem != query_props->allowed_list.front()) {
                element_list_text.append(", ");
            }
            auto symbol = atomic_number_to_symbol(static_cast<int>(cur_elem));
            element_list_text.append(symbol);
        }
        if (element_list_text.empty()) {
            element_list_text = "Empty";
        }

        os << "AtomQueryProperties:\n"
           << "\tQuery type: " << static_cast<int>(query_props->query_type)
           << "\n"
           << "\tElement: "
           << atomic_number_to_symbol(static_cast<int>(query_props->element))
           << "\n"
           << "\tElement list: " << element_list_text << "\n"
           << "\tIsotope: " << output_optional(query_props->isotope) << "\n"
           << "\tCharge: " << output_optional(query_props->charge) << "\n"
           << "\tUnpaired electrons: "
           << output_optional(query_props->unpaired_electrons) << "\n"
           << "\tWildcard: " << static_cast<int>(query_props->wildcard) << "\n"
           << "\tR-group: " << static_cast<int>(query_props->r_group) << "\n"
           << "\tSMARTS query: \"" << query_props->smarts_query << "\"\n"
           << "\tTotal H type: " << static_cast<int>(query_props->total_h_type)
           << "\n"
           << "\tTotal H exact value: " << query_props->total_h_exact_val
           << "\n"
           << "\tNum connections: "
           << output_optional(query_props->num_connections) << "\n"
           << "\tAromaticity: " << static_cast<int>(query_props->aromaticity)
           << "\n"
           << "\tRing count type: "
           << static_cast<int>(query_props->ring_count_type) << "\n"
           << "\tRing count exact value: " << query_props->ring_count_exact_val
           << "\n"
           << "\tRing bond count type: "
           << static_cast<int>(query_props->ring_bond_count_type) << "\n"
           << "\tRing bond count exact value: "
           << query_props->ring_bond_count_exact_val << "\n"
           << "\tSmallest ring size: "
           << output_optional(query_props->smallest_ring_size) << "\n";
    }
    return os;
}

bool AtomQueryProperties::hasPropertiesBeyondQueryType() const
{
    // General properties excluding enhanced_stereo (stereo doesn't affect
    // whether we can display allowed lists with element symbols - SKETCH-2487)
    const AtomQueryPropertyList GENERAL_PROPS_EXCL_STEREO = {
        &AtomQueryProperties::isotope,
        &AtomQueryProperties::charge,
        &AtomQueryProperties::unpaired_electrons,
    };
    auto default_props = AtomQueryProperties();
    return !compare_query_properties(
        this, &default_props,
        concat_query_props(
            {GENERAL_PROPS_EXCL_STEREO, ADVANCED_QUERY_PROPERTIES}));
}

bool AtomQueryProperties::hasAdvancedProperties() const
{
    auto default_props = AtomQueryProperties();
    return !compare_query_properties(this, &default_props,
                                     ADVANCED_QUERY_PROPERTIES);
}

/**
 * Throw an UnrecognizedQueryError if the given query is negated
 */
static void throw_if_negated(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    if (query->getNegation()) {
        throw UnrecognizedQueryError(
            "Cannot negate a query with this description");
    }
}

/**
 * If the given query is an equality query (i.e. <atom property> equals <n>, as
 * opposed to a less than or greater than query), then return the matching
 * value.
 * @throw UnrecognizedQueryError if the given query is not an equality query
 */
static int
get_value_for_equality_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    if (auto* equality_query =
            dynamic_cast<const RDKit::ATOM_EQUALS_QUERY*>(query)) {
        return equality_query->getVal();
    } else {
        throw UnrecognizedQueryError(
            "Can only parse equality queries for this property");
    }
}

/**
 * @overload A version of this function that handles optional property values
 */
template <typename T> static void
check_value_for_conflicts(T new_val, std::optional<T> existing_val,
                          std::string desc,
                          std::unordered_set<std::string>& seen_descriptions)
{
    if (existing_val.has_value()) {
        check_value_for_conflicts<T>(new_val, *existing_val, desc,
                                     seen_descriptions);
    } else {
        seen_descriptions.insert(desc);
    }
}

/**
 * When we read in new portion of a query, make sure that it doesn't conflict
 * with any of the other queries on the same atom (e.g. make sure that the
 * query doesn't specify anything like "this atom is in exactly one ring AND
 * this atom is in exactly two rings") and record that we've seen the new query.
 * @param new_val The query value that was just read in
 * @param existing_val The query value for this property that was previously
 * read in, if any. This value is only used if `desc` appears in
 * `seen_descriptions`.
 * @param desc The description of the current query. This is typically the
 * return value from QUERYATOM_QUERY::getDescription(). However, if multiple
 * query types can affect the the same property (e.g. AtomIsAliphatic and
 * AtomIsAromatic both affect AtomQueryProperties::aromaticity), then use one of
 * the constants at the top of this file (e.g. ATOM_AROMATICITY) instead.
 * @param seen_descriptions The set of all `desc` values that we've already seen
 * for this atom. As part of this function, `desc` will be added to this set.
 */
template <typename T> static void
check_value_for_conflicts(T new_val, T existing_val, std::string desc,
                          std::unordered_set<std::string>& seen_descriptions)
{
    if (seen_descriptions.count(desc)) {
        if (new_val != existing_val) {
            throw UnrecognizedQueryError("Query specifies conflicting values");
        }
    } else {
        seen_descriptions.insert(desc);
    }
}

/**
 * Convert a signed integer to unsigned
 * @throw UnrecognizedQueryError if the value is negative
 */
static unsigned int as_unsigned(int val)
{
    if (val < 0) {
        throw UnrecognizedQueryError(
            "Negative values not allowed for this query type");
    }
    return static_cast<unsigned int>(val);
}

/**
 * If the given query is a wildcard query, add it to query_props. If not, do
 * nothing.
 * @param query A query instance that might be a wildcard query
 * @param query_props The query properties that we have already read in. If
 * `query` is a wildcard query, then `query_type` and `wildcard` will be set on
 * `query_props`.
 * @param seen_descriptions The set of all query description values that we've
 * already seen for this atom. As part of this function, `QUERY_TYPE` will be
 * added to this set if and only if this is a wildcard query.
 * @return whether the query was a wildcard query
 */
static bool process_possible_wildcard_query(
    const RDKit::Atom::QUERYATOM_QUERY* const query,
    std::shared_ptr<AtomQueryProperties> query_props,
    std::unordered_set<std::string>& seen_descriptions)
{
    auto type_label = query->getTypeLabel();
    if (!TYPE_LABEL_TO_ATOM_QUERY.count(type_label)) {
        return false;
    }
    auto wildcard = TYPE_LABEL_TO_ATOM_QUERY.at(type_label);
    auto query_type = QueryType::WILDCARD;
    std::pair<QueryType, AtomQuery> new_vals = {query_type, wildcard};
    std::pair<QueryType, AtomQuery> existing_vals = {query_props->query_type,
                                                     query_props->wildcard};
    check_value_for_conflicts(new_vals, existing_vals, QUERY_TYPE,
                              seen_descriptions);
    query_props->query_type = query_type;
    query_props->wildcard = wildcard;
    return true;
}

/**
 * @return the element for the specified atomic number
 * @throw UnrecognizedQueryError if the atomic number is invalid
 */
static Element get_element_for_atomic_num(int atomic_num)
{
    if (!is_atomic_number(atomic_num)) {
        throw UnrecognizedQueryError("Invalid atomic number");
    }
    return Element(atomic_num);
}

/**
 * Read in the element and aromaticity information from an AtomType query.
 * Aromatic queries are indicated by a query value of 1000 plus the atomic
 * number.  Values of less than 1000 indicate an aliphatic query.
 */
static std::tuple<Element, QueryAromaticity>
get_element_and_aromaticity_from_atom_type_query(
    const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    auto atom_type_num = get_value_for_equality_query(query);
    auto aromaticity = atom_type_num > ATOM_TYPE_AROMATIC_OFFSET
                           ? QueryAromaticity::AROMATIC
                           : QueryAromaticity::ALIPHATIC;
    auto atomic_number = atom_type_num % ATOM_TYPE_AROMATIC_OFFSET;
    auto element = get_element_for_atomic_num(atomic_number);
    return {element, aromaticity};
}

/**
 * Parse all descendents of an AtomOr query to determine the allowed elements
 * and the aromaticity specified, if any.
 * @return A tuple of
 *   - the aromaticity, if one is specified, or QueryAromaticity::ANY if not
 *   - the vector of specified elements, returned in the order they appeared in
 *     the query.  The elements of this vector are guaranteed to be unique (i.e.
 *     no duplicates)
 *   - the set of specified elements.  This set contains the same elements as
 *     the above vector.
 * @throw UnrecognizedQueryError if the query specifies anything other than
 * allowed elements or aromaticity, or if the aromaticity differs for different
 * elements
 */
static std::tuple<QueryAromaticity, std::vector<Element>,
                  std::unordered_set<Element>>
parse_children_of_or_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    auto aromaticity = QueryAromaticity::ANY;
    bool seen_first_element = false;
    std::vector<Element> elements;
    std::unordered_set<Element> elements_set;
    std::deque<std::shared_ptr<RDKit::Atom::QUERYATOM_QUERY>> queries;
    auto append_children_to_queries =
        [&queries](const RDKit::Atom::QUERYATOM_QUERY* const query) {
            // in order for the element list to match the input order, we have
            // to prepend the children, in order, to our deque
            std::vector<std::shared_ptr<RDKit::Atom::QUERYATOM_QUERY>> children(
                query->beginChildren(), query->endChildren());
            for (auto child_it = children.rbegin(); child_it != children.rend();
                 ++child_it) {
                queries.push_front(*child_it);
            }
        };
    // add the element to the elements vector if it's not already on the list
    auto record_element = [&elements, &elements_set](Element cur_element) {
        if (!elements_set.count(cur_element)) {
            elements.push_back(cur_element);
            elements_set.insert(cur_element);
        }
    };

    append_children_to_queries(query);
    while (!queries.empty()) {
        auto cur_child = queries.front();
        queries.pop_front();
        throw_if_negated(cur_child.get());
        auto desc = cur_child->getDescription();
        if (desc == "AtomAtomicNum") {
            seen_first_element = true;
            if (aromaticity != QueryAromaticity::ANY) {
                // throw an exception if the AtomOr query mixes AtomAtomicNum
                // children (which can't specify an aromaticity) with AtomType
                // children (which always specify an aromaticity)
                throw UnrecognizedQueryError(
                    "Cannot specify different aromaticity for different "
                    "elements");
            }
            auto atomic_num = get_value_for_equality_query(cur_child.get());
            auto cur_element = get_element_for_atomic_num(atomic_num);
            record_element(cur_element);
        } else if (desc == "AtomType") {
            auto [cur_element, cur_aromaticity] =
                get_element_and_aromaticity_from_atom_type_query(
                    cur_child.get());
            if (!seen_first_element) {
                aromaticity = cur_aromaticity;
                seen_first_element = true;
            } else if (cur_aromaticity != aromaticity) {
                throw UnrecognizedQueryError(
                    "Cannot specify different aromaticity for different "
                    "elements");
            }
            record_element(cur_element);
        } else if (desc == "AtomOr") {
            append_children_to_queries(cur_child.get());
        } else {
            throw UnrecognizedQueryError(
                "Cannot parse AtomOr queries that specify properties other "
                "than elements");
        }
    }
    return {aromaticity, elements, elements_set};
}

/**
 * Store a query for the specified element in the AtomQueryProperties object
 * @param[in] element The queried element
 * @param[in] is_negated Whether the query was negated
 * @param[in,out] query_props The object to store the query in
 * @param[in,out] seen_descriptions The set of all query description values that
 * we've already seen for this atom. This set will be updated to note this
 * element query
 * @throw UnrecognizedQueryError if this query conflicts with a previous element
 * query (i.e. if this atom's query describes both a required specified element
 * *and* a disallowed element list)
 */
static void
process_element_query(const Element element, const bool is_negated,
                      std::shared_ptr<AtomQueryProperties> query_props,
                      std::unordered_set<std::string>& seen_descriptions)
{
    if (!is_negated) {
        auto query_type = QueryType::SPECIFIC_ELEMENT;
        std::pair<QueryType, Element> new_vals = {query_type, element};
        std::pair<QueryType, Element> existing_vals = {query_props->query_type,
                                                       query_props->element};
        check_value_for_conflicts(new_vals, existing_vals, QUERY_TYPE,
                                  seen_descriptions);
        query_props->query_type = query_type;
        query_props->element = element;
    } else {
        auto query_type = QueryType::NOT_ALLOWED_LIST;
        check_value_for_conflicts(query_type, query_props->query_type,
                                  QUERY_TYPE, seen_descriptions);
        query_props->query_type = query_type;
        query_props->allowed_list.push_back(element);
    }
    // an ALLOWED_LIST query will always be a child of an AtomOr query, so
    // that's handled in parse_children_of_or_query above
}

/**
 * Read in all of the properties of the given query, recursing into children
 * queries if any are present. Any properties that are read in will be added to
 * query_props.
 *
 * This function should only be called from itself or from read_query. Other
 * callers should use read_query instead.
 *
 * @param[in] query The query to read in
 * @param[in,out] query_props The query properties instance to update with the
 * read-in properties
 * @param[in,out] seen_descriptions The set of all query description values that
 * we've already seen for this atom. This set will be updated as we read in the
 * given query.
 */
static void
read_query_recursive(const RDKit::Atom::QUERYATOM_QUERY* const query,
                     std::shared_ptr<AtomQueryProperties> query_props,
                     std::unordered_set<std::string>& seen_descriptions)
{
    auto desc = query->getDescription();
    if (desc == "AtomAnd") {
        for (auto it = query->beginChildren(); it != query->endChildren();
             ++it) {
            read_query_recursive(it->get(), query_props, seen_descriptions);
        }
    } else if (desc == "AtomOr") {
        bool is_wildcard = process_possible_wildcard_query(query, query_props,
                                                           seen_descriptions);
        if (!is_wildcard) {
            auto query_type = query->getNegation() ? QueryType::NOT_ALLOWED_LIST
                                                   : QueryType::ALLOWED_LIST;
            auto [aromaticity, elements, elements_set] =
                parse_children_of_or_query(query);

            // as long as this query isn't empty, store the query type and the
            // allowed or disallowed list
            if (!elements.empty()) {
                if (query_type == QueryType::ALLOWED_LIST) {
                    std::pair<QueryType, std::unordered_set<Element>> new_vals =
                        {query_type, elements_set};
                    std::unordered_set existing_elements_set(
                        query_props->allowed_list.begin(),
                        query_props->allowed_list.end());
                    std::pair<QueryType, std::unordered_set<Element>>
                        existing_vals = {query_props->query_type,
                                         existing_elements_set};
                    // If we have two allowed lists in the query, then both
                    // lists have to specify the exact same elements, since
                    // there's no way to match a query like "#6,#7;+2;#8,#9".
                    // An element can't match both "carbon or nitrogen" and
                    // "oxygen or fluorine" at the same time.
                    check_value_for_conflicts(new_vals, existing_vals,
                                              QUERY_TYPE, seen_descriptions);
                    query_props->query_type = query_type;
                    query_props->allowed_list = elements;
                } else {
                    // If we have two disallowed lists in the query, then the
                    // lists *don't* need to specify the same elements; we can
                    // just append the two lists together.  In a query like
                    // "!#6&!#7&+2&!#8&!#9", then it doesn't matter that the +2
                    // charge term breaks up the "not carbon and not nitrogen"
                    // list from the "not oxygen and not fluorine" list; any
                    // element other than C, N, O, or F will match both
                    // disallowed lists.
                    check_value_for_conflicts(query_type,
                                              query_props->query_type,
                                              QUERY_TYPE, seen_descriptions);
                    query_props->query_type = query_type;
                    std::copy(elements.begin(), elements.end(),
                              query_props->allowed_list.end());
                }
            }
            // if the child queries were all AtomType queries (which specify
            // both atomic number and aromaticity), store the aromaticity
            if (aromaticity != QueryAromaticity::ANY) {
                check_value_for_conflicts(aromaticity, query_props->aromaticity,
                                          ATOM_AROMATICITY, seen_descriptions);
                query_props->aromaticity = aromaticity;
            }
        }
    } else if (desc == "AtomAtomicNum") {
        // some wildcards are stored as negated atomic number queries (e.g. A is
        // not H, QH is not C)
        bool is_wildcard = process_possible_wildcard_query(query, query_props,
                                                           seen_descriptions);
        if (!is_wildcard) {
            auto atomic_number = get_value_for_equality_query(query);
            if (!is_atomic_number(atomic_number)) {
                throw UnrecognizedQueryError("Invalid atomic number");
            }
            auto element = Element(atomic_number);
            auto is_negated = query->getNegation();
            process_element_query(element, is_negated, query_props,
                                  seen_descriptions);
        }
    } else if (desc == "AtomFormalCharge") {
        throw_if_negated(query);
        auto val = get_value_for_equality_query(query);
        check_value_for_conflicts(val, query_props->charge, desc,
                                  seen_descriptions);
        query_props->charge = val;
    } else if (desc == "AtomHCount") {
        auto val = as_unsigned(get_value_for_equality_query(query));
        auto new_type = QueryCount::EXACTLY;
        if (query->getNegation()) {
            if (val != 0) {
                // this query is checking whether the number of hydrogens is not
                // equal to something other than zero, which can't be
                // represented in the atom properties dialog
                throw UnrecognizedQueryError("AtomHCount queries may only be "
                                             "an exact value or not zero.");
            }
            new_type = QueryCount::POSITIVE;
        }
        std::pair<QueryCount, unsigned int> new_total_h_vals = {new_type, val};
        std::pair<QueryCount, unsigned int> existing_total_h_vals = {
            query_props->total_h_type, query_props->total_h_exact_val};
        check_value_for_conflicts(new_total_h_vals, existing_total_h_vals, desc,
                                  seen_descriptions);
        query_props->total_h_type = new_type;
        query_props->total_h_exact_val = val;
    } else if (desc == "AtomInNRings") {
        throw_if_negated(query);
        auto new_type = QueryCount::EXACTLY;
        unsigned int new_exact_val = get_value_for_equality_query(query);
        std::pair<QueryCount, unsigned int> new_ring_vals = {new_type,
                                                             new_exact_val};
        std::pair<QueryCount, unsigned int> existing_ring_vals = {
            query_props->ring_count_type, query_props->ring_count_exact_val};
        check_value_for_conflicts(new_ring_vals, existing_ring_vals,
                                  ATOM_RING_COUNT, seen_descriptions);
        query_props->ring_count_type = new_type;
        query_props->ring_count_exact_val = new_exact_val;
    } else if (desc == "AtomInRing") {
        auto new_type =
            query->getNegation() ? QueryCount::EXACTLY : QueryCount::POSITIVE;
        unsigned int new_exact_val = 0;
        std::pair<QueryCount, unsigned int> new_ring_vals = {new_type,
                                                             new_exact_val};
        std::pair<QueryCount, unsigned int> existing_ring_vals = {
            query_props->ring_count_type, query_props->ring_count_exact_val};
        check_value_for_conflicts(new_ring_vals, existing_ring_vals,
                                  ATOM_RING_COUNT, seen_descriptions);
        query_props->ring_count_type = new_type;
        query_props->ring_count_exact_val = new_exact_val;
    } else if (desc == "AtomIsAliphatic") {
        auto val = query->getNegation() ? QueryAromaticity::AROMATIC
                                        : QueryAromaticity::ALIPHATIC;
        check_value_for_conflicts(val, query_props->aromaticity,
                                  ATOM_AROMATICITY, seen_descriptions);
        query_props->aromaticity = val;
    } else if (desc == "AtomIsAromatic") {
        auto val = query->getNegation() ? QueryAromaticity::ALIPHATIC
                                        : QueryAromaticity::AROMATIC;
        check_value_for_conflicts(val, query_props->aromaticity,
                                  ATOM_AROMATICITY, seen_descriptions);
        query_props->aromaticity = val;
    } else if (desc == "AtomIsotope") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->isotope, desc,
                                  seen_descriptions);
        query_props->isotope = val;
    } else if (desc == "AtomMinRingSize") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->smallest_ring_size, desc,
                                  seen_descriptions);
        query_props->smallest_ring_size = val;
    } else if (desc == "AtomNull") {
        // AH wildcard queries are stored as a null query
        process_possible_wildcard_query(query, query_props, seen_descriptions);
        // this could also be an R-group, but that's handled in
        // read_properties_from_atom since the R-group number is stored in an
        // atom property, not in the query
    } else if (desc == "AtomNumRadicalElectrons") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->unpaired_electrons, desc,
                                  seen_descriptions);
        query_props->unpaired_electrons = val;
    } else if (desc == "AtomRingBondCount") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->ring_bond_count_exact_val,
                                  desc, seen_descriptions);
        query_props->ring_bond_count_type = QueryCount::EXACTLY;
        query_props->ring_bond_count_exact_val = val;
    } else if (desc == "AtomTotalDegree") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->num_connections, desc,
                                  seen_descriptions);
        query_props->num_connections = val;
    } else if (desc == "AtomType") {
        bool is_negated = query->getNegation();
        auto [element, aromaticity] =
            get_element_and_aromaticity_from_atom_type_query(query);
        process_element_query(element, is_negated, query_props,
                              seen_descriptions);
        check_value_for_conflicts(aromaticity, query_props->aromaticity,
                                  ATOM_AROMATICITY, seen_descriptions);
        query_props->aromaticity = aromaticity;
    } else {
        throw UnrecognizedQueryError("Description not recognized");
    }
}

/**
 * Convert the given query into a SMARTS string
 */
std::string
get_smarts_for_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    auto* query_atom = new RDKit::QueryAtom();
    query_atom->setQuery(query->copy());
    auto smarts = RDKit::SmartsWrite::GetAtomSmarts(query_atom);
    delete query_atom;
    return smarts;
}

/**
 * Read in all of the properties of the given query.
 */
std::shared_ptr<AtomQueryProperties>
read_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    auto query_props = std::make_shared<AtomQueryProperties>();
    std::unordered_set<std::string> seen_descriptions;
    try {
        read_query_recursive(query, query_props, seen_descriptions);

        // if we haven't read any query criteria, then the query represents the
        // AH wildcard.  (We should only hit this with a query that was read in
        // from SMARTS.  Otherwise, the query would have a type label that
        // would've been recognized by read_query_recursive.)
        if (seen_descriptions.empty()) {
            query_props->query_type = QueryType::WILDCARD;
            query_props->wildcard = AtomQuery::AH;
        }
    } catch (const UnrecognizedQueryError&) {
        // we couldn't parse the query into something that could be represented
        // in the dialog's property field, so report it as a SMARTS query
        // instead
        auto smarts = get_smarts_for_query(query);
        // reset any properties that we did manage to read
        query_props.reset(new AtomQueryProperties());

        query_props->query_type = QueryType::SMARTS;
        query_props->smarts_query = smarts;
    }
    return query_props;
}

/**
 * Read in all properties from the given atom
 */
std::shared_ptr<AbstractAtomProperties>
read_properties_from_atom(const RDKit::Atom* const atom)
{
    std::shared_ptr<AbstractAtomProperties> props = nullptr;
    auto r_group_num = rdkit_extensions::get_r_group_number(atom);
    if (r_group_num.has_value()) {
        // this is an R-group, which is stored as a non-query Atom, but
        // represented in the dialog as a query
        auto* query_props = new AtomQueryProperties();
        props.reset(query_props);
        query_props->query_type = QueryType::RGROUP;
        query_props->r_group = *r_group_num;
    } else if (!atom->hasQuery()) {
        auto* atom_props = new AtomProperties();
        props.reset(atom_props);
        atom_props->element = Element(atom->getAtomicNum());
        auto isotope = atom->getIsotope();
        if (isotope == 0) {
            atom_props->isotope = std::nullopt;
        } else {
            atom_props->isotope = isotope;
        }
        atom_props->charge = atom->getFormalCharge();
        atom_props->unpaired_electrons = atom->getNumRadicalElectrons();
    } else {
        auto* query_atom = static_cast<const RDKit::QueryAtom*>(atom);
        props = read_query(query_atom->getQuery());
    }
    props->enhanced_stereo = get_enhanced_stereo_for_atom(atom);
    return props;
}

/**
 * Add the query to the atom by calling expandQuery (if the atom already has a
 * query) or setQuery (if it doesn't)
 * @param[in] query The query to add
 * @param[in,out] query_atom The atom to add the query to
 * @param[in] logical_op The logical operator (i.e. AND or OR) to use if calling
 * expandQuery.  Note that this argument will be ignored if the atom doesn't yet
 * have a query set on it.
 */
static void add_query_to_atom(
    RDKit::Atom::QUERYATOM_QUERY* const query,
    const std::shared_ptr<RDKit::Atom> query_atom,
    const Queries::CompositeQueryType logical_op = Queries::COMPOSITE_AND)
{
    if (!query_atom->hasQuery()) {
        query_atom->setQuery(query);
    } else {
        query_atom->expandQuery(query, logical_op);
    }
}

QString get_label_from_atom_query(AtomQuery query)
{
    for (auto [label, atom_query] : TYPE_LABEL_TO_ATOM_QUERY) {
        if (atom_query == query) {
            return QString::fromStdString(label);
        }
    }
    return "*";
}

QString get_atomic_symbol_from_element(Element element)
{
    auto symbol = atomic_number_to_symbol(static_cast<unsigned int>(element));
    return QString::fromStdString(symbol);
}

QString join_all_atomic_symbols(std::vector<Element> elements,
                                QString separator)
{
    QStringList atomic_symbols;
    std::transform(elements.begin(), elements.end(),
                   std::back_inserter(atomic_symbols),
                   get_atomic_symbol_from_element);
    return atomic_symbols.join(separator);
}

/**
 * Create an atom instance using the given query properties. The newly created
 * atom will have any queries and/or atom properties necessary to describe the
 * query_props->query_type value.
 * @param query_props The query properties to use
 * @return The newly created atom.  Note that this will be a QueryAtom instance
 * *unless* the query_props specified an R-group query. In that case, it will be
 * a non-query Atom instance with the appropriate R-group atom label.
 */
static std::shared_ptr<RDKit::Atom> create_query_atom_for_query_type_property(
    const AtomQueryProperties* const query_props)
{
    std::shared_ptr<RDKit::Atom> query_atom;
    // first, create the atom
    if (query_props->query_type == QueryType::RGROUP) {
        query_atom = rdkit_extensions::make_new_r_group(query_props->r_group);
    } else {
        query_atom = std::make_shared<RDKit::QueryAtom>();
    }

    bool allowed = false;
    switch (query_props->query_type) {
        case QueryType::ALLOWED_LIST:
            if (query_props->allowed_list.size() == 1) {
                auto atomic_num =
                    static_cast<int>(*query_props->allowed_list.begin());
                query_atom->setAtomicNum(atomic_num);
            }
            allowed = true;
            [[fallthrough]];
        case QueryType::NOT_ALLOWED_LIST: {
            auto logical_op =
                allowed ? Queries::COMPOSITE_OR : Queries::COMPOSITE_AND;
            for (auto element : query_props->allowed_list) {
                auto* cur_query =
                    RDKit::makeAtomNumQuery(static_cast<int>(element));
                cur_query->setNegation(!allowed);
                add_query_to_atom(cur_query, query_atom, logical_op);
            }
            break;
        }
        case QueryType::WILDCARD: {
            auto wildcard_query =
                ATOM_TOOL_QUERY_MAP.at(query_props->wildcard)();
            add_query_to_atom(wildcard_query, query_atom);
            break;
        }
        case QueryType::SPECIFIC_ELEMENT: {
            auto atomic_num = static_cast<int>(query_props->element);
            query_atom->setAtomicNum(atomic_num);
            add_query_to_atom(RDKit::makeAtomNumQuery(atomic_num), query_atom);
            break;
        }
        case QueryType::RGROUP: {
            // nothing to do here, since we already set the R-group when
            // creating the atom
            break;
        }
        case QueryType::SMARTS: {
            auto smarts = boost::trim_copy(query_props->smarts_query);
            auto* smarts_atom = RDKit::SmartsToAtom(smarts);
            if (smarts_atom == nullptr) {
                // we should never hit this since the dialog disables the OK
                // button when the SMARTS is invalid
                throw std::runtime_error("Invalid SMARTS pattern: " + smarts);
            }
            auto* query = smarts_atom->getQuery()->copy();
            // if the atomic number isn't right, RDKit may ignore the query?
            query_atom->setAtomicNum(smarts_atom->getAtomicNum());
            delete smarts_atom;
            add_query_to_atom(query, query_atom);
        }
    }
    return query_atom;
}

/**
 * Update the atom with all properties that were specified in the General tab.
 * @param[in] query_props The query properties to transfer to query_atom
 * @param[in,out] query_atom The query atom to update
 */
static void
update_atom_for_general_properties(const AtomQueryProperties* const query_props,
                                   std::shared_ptr<RDKit::Atom> query_atom)
{
    if (query_props->isotope.has_value()) {
        add_query_to_atom(RDKit::makeAtomIsotopeQuery(*query_props->isotope),
                          query_atom);
        query_atom->setIsotope(*query_props->isotope);
    }
    if (query_props->charge.has_value()) {
        auto query = RDKit::makeAtomFormalChargeQuery(*query_props->charge);
        add_query_to_atom(query, query_atom);
        query_atom->setFormalCharge(*query_props->charge);
    }
    if (query_props->unpaired_electrons.has_value()) {
        auto* query = RDKit::makeAtomNumRadicalElectronsQuery(
            *query_props->unpaired_electrons);
        add_query_to_atom(query, query_atom);
        query_atom->setNumRadicalElectrons(*query_props->unpaired_electrons);
    }
}

/**
 * Update the atom with all properties that were specified in the Advanced tab.
 * @param[in] query_props The query properties to transfer to query_atom
 * @param[in,out] query_atom The query atom to update
 */
static void update_atom_for_advanced_properties(
    const AtomQueryProperties* const query_props,
    std::shared_ptr<RDKit::Atom> query_atom)
{
    if (query_props->total_h_type != QueryCount::ANY) {
        int query_val = 0;
        if (query_props->total_h_type == QueryCount::EXACTLY) {
            query_val = query_props->total_h_exact_val;
            query_atom->setNumExplicitHs(query_props->total_h_exact_val);
        }
        auto* query = RDKit::makeAtomHCountQuery(query_val);
        query->setNegation(query_props->total_h_type == QueryCount::POSITIVE);
        add_query_to_atom(query, query_atom);
    }

    if (query_props->num_connections.has_value()) {
        // the BIOVIA specification for NCNN for exactly 0 connections is -1
        // (0 means unspecified), so we need to account for it when setting the
        // molNumBonds property
        int subst_count = *query_props->num_connections == 0
                              ? -1
                              : *query_props->num_connections;
        query_atom->setProp(RDKit::common_properties::molSubstCount,
                            subst_count);
        auto* query =
            RDKit::makeAtomTotalDegreeQuery(*query_props->num_connections);
        add_query_to_atom(query, query_atom);
    }

    if (query_props->aromaticity == QueryAromaticity::AROMATIC) {
        auto* query = RDKit::makeAtomAromaticQuery();
        add_query_to_atom(query, query_atom);
        query_atom->setIsAromatic(true);
    } else if (query_props->aromaticity == QueryAromaticity::ALIPHATIC) {
        auto* query = RDKit::makeAtomAliphaticQuery();
        add_query_to_atom(query, query_atom);
        query_atom->setIsAromatic(false);
    }

    if (query_props->ring_count_type == QueryCount::POSITIVE) {
        auto* query = RDKit::makeAtomInRingQuery();
        add_query_to_atom(query, query_atom);
    } else if (query_props->ring_count_type == QueryCount::EXACTLY) {
        auto* query =
            RDKit::makeAtomInNRingsQuery(query_props->ring_count_exact_val);
        add_query_to_atom(query, query_atom);
    }

    if (query_props->ring_bond_count_type == QueryCount::EXACTLY) {
        //  the BIOVIA specification for RBCNT for exactly 0 ring bonds is -1 (0
        //  means unspecified), so we need to account for it when setting the
        //  molRingBondCount property
        int rbcnt = query_props->ring_bond_count_exact_val == 0
                        ? -1
                        : query_props->ring_bond_count_exact_val;
        query_atom->setProp(RDKit::common_properties::molRingBondCount, rbcnt);
        auto* query = RDKit::makeAtomRingBondCountQuery(
            query_props->ring_bond_count_exact_val);
        add_query_to_atom(query, query_atom);
    }
    if (query_props->smallest_ring_size.has_value()) {
        auto* query =
            RDKit::makeAtomMinRingSizeQuery(*query_props->smallest_ring_size);
        add_query_to_atom(query, query_atom);
    }
}

// TODO: according to ChemicalKnowledge.cpp, product atoms can't have
//       allowed/disallowed list element queries, so it silently discards them.
//       (See SKETCH-843 and ENUM-402.)  Should probably handle this at the
//       export level.
/**
 * Create an RDKit atom with the given properties.
 */
std::pair<std::shared_ptr<RDKit::Atom>, std::optional<EnhancedStereo>>
create_atom_with_properties(
    const std::shared_ptr<AbstractAtomProperties> properties)
{
    std::shared_ptr<RDKit::Atom> atom = nullptr;
    if (!properties->isQuery()) {
        auto* atom_props = static_cast<const AtomProperties*>(properties.get());
        atom = std::make_shared<RDKit::Atom>(
            static_cast<int>(atom_props->element));
        atom->setIsotope(atom_props->isotope.value_or(0));
        atom->setFormalCharge(atom_props->charge);
        atom->setNumRadicalElectrons(atom_props->unpaired_electrons);
    } else {
        auto* query_props =
            static_cast<const AtomQueryProperties*>(properties.get());
        atom = create_query_atom_for_query_type_property(query_props);
        update_atom_for_general_properties(query_props, atom);
        update_atom_for_advanced_properties(query_props, atom);
    }
    return {atom, properties->enhanced_stereo};
}

std::shared_ptr<AtomQueryProperties>
get_only_advanced_properties(std::shared_ptr<AtomQueryProperties> property)
{
    auto advanced_props =
        get_subset_of_properties(property.get(), ADVANCED_QUERY_PROPERTIES);
    return std::make_shared<AtomQueryProperties>(advanced_props);
}
} // namespace schrodinger::sketcher
