#include "schrodinger/sketcher/rdkit/atom_properties.h"

#include <unordered_set>

#include <boost/algorithm/string.hpp>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
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

namespace
{

class UnrecognizedQueryError : public std::runtime_error
{
  public:
    using std::runtime_error::runtime_error;
};

} // namespace

bool EnhancedStereo::operator==(const EnhancedStereo& other) const
{
    return (this->type == other.type && this->group_id == other.group_id);
}

bool EnhancedStereo::operator!=(const EnhancedStereo& other) const
{
    return !(*this == other);
}

bool AbstractAtomProperties::operator==(
    const AbstractAtomProperties& other) const
{
    if (isQuery() != other.isQuery() || this->element != other.element ||
        this->isotope != other.isotope ||
        this->enhanced_stereo != other.enhanced_stereo) {
        return false;
    } else if (isQuery()) {
        auto* this_query = static_cast<const AtomQueryProperties*>(this);
        auto* other_query = static_cast<const AtomQueryProperties*>(&other);
        return this_query->charge == other_query->charge &&
               this_query->unpaired_electrons ==
                   other_query->unpaired_electrons &&
               this_query->query_type == other_query->query_type &&
               this_query->allowed_list == other_query->allowed_list &&
               this_query->wildcard == other_query->wildcard &&
               this_query->r_group == other_query->r_group &&
               this_query->total_h == other_query->total_h &&
               this_query->num_connections == other_query->num_connections &&
               this_query->aromaticity == other_query->aromaticity &&
               this_query->ring_count_type == other_query->ring_count_type &&
               this_query->ring_count_exact_val ==
                   other_query->ring_count_exact_val &&
               this_query->ring_bond_count_type ==
                   other_query->ring_bond_count_type &&
               this_query->ring_bond_count_exact_val ==
                   other_query->ring_bond_count_exact_val &&
               this_query->smallest_ring_size ==
                   other_query->smallest_ring_size &&
               this_query->smarts_query == other_query->smarts_query;
    } else { // not a query
        auto* this_atom = static_cast<const AtomProperties*>(this);
        auto* other_atom = static_cast<const AtomProperties*>(&other);
        return this_atom->charge == other_atom->charge &&
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
           << atomic_number_to_name(static_cast<int>(atom_props->element))
           << "\n"
           << "\tIsotope: " << output_optional(atom_props->isotope) << "\n"
           << "\tCharge: " << atom_props->charge << "\n"
           << "\tUnpaired electrons: " << atom_props->unpaired_electrons
           << "\n";
    } else {
        auto* query_props = static_cast<const AtomQueryProperties*>(&props);
        os << "AtomQueryProperties:\n"
           << "\tQuery Type: " << static_cast<int>(query_props->query_type)
           << "\n"
           // TODO: allowed_list
           << "\tElement: "
           << atomic_number_to_name(static_cast<int>(query_props->element))
           << "\n"
           << "\tIsotope: " << output_optional(query_props->isotope) << "\n"
           << "\tCharge: " << output_optional(query_props->charge) << "\n"
           << "\tUnpaired electrons: "
           << output_optional(query_props->unpaired_electrons) << "\n"
           << "\tWildcard: " << static_cast<int>(query_props->wildcard) << "\n"
           << "\tR-group: " << static_cast<int>(query_props->r_group) << "\n"
           << "\tTotal H: " << output_optional(query_props->total_h) << "\n"
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

/**
 * Read the enhanced stereo data, if any, for an atom.
 */
static std::optional<EnhancedStereo>
read_enhanced_stereo_from_atom(const RDKit::Atom* const atom)
{
    auto& mol = atom->getOwningMol();
    for (const auto& stereo_group : mol.getStereoGroups()) {
        const auto& group_atoms = stereo_group.getAtoms();
        if (std::find(group_atoms.begin(), group_atoms.end(), atom) !=
            group_atoms.end()) {
            EnhancedStereo enhanced_stereo;
            enhanced_stereo.type = stereo_group.getGroupType();
            enhanced_stereo.group_id = stereo_group.getWriteId();
            return enhanced_stereo;
        }
    }
    // this atom doesn't appear in any stereo groups
    return std::nullopt;
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

/*!
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

static std::pair<std::vector<Element>, QueryAromaticity>
parse_children_of_or_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    // TODO: only allow children to specify AtomType or AtomAtomicNum, or
    //       AtomAnd combinations of these (so long as the AtomType is
    //       consistent)
    return {};
}

/**
 * Read in all of the properties of the given query, recursing into children
 * queries if any are present. Any properties that are read in will be added to
 * query_props.
 *
 * This function should only be called from itself or from read_query. Other
 * callers should use read_query instead.
 *
 * @param query The query to read in
 * @param query_props The query properties instance to update with the read-in
 * properties
 * @param seen_descriptions The set of all query description values that we've
 * already seen for this atom. This set will be updated as we read in the given
 * query.
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
            auto [elements, aromaticity] = parse_children_of_or_query(query);
            // auto query_type = query->getNegation() ?
            // QueryType::NOT_ALLOWED_LIST
            //                                        : QueryType::ALLOWED_LIST;
            // TODO: convert element lists to sets so we don't compare order
            // TODO: store query type and elements

            if (aromaticity != QueryAromaticity::ANY) {
                check_value_for_conflicts(aromaticity, query_props->aromaticity,
                                          ATOM_AROMATICITY, seen_descriptions);
                query_props->aromaticity = aromaticity;
            }
        }
    } else if (desc == "AtomAtomicNum") {
        bool is_wildcard = process_possible_wildcard_query(query, query_props,
                                                           seen_descriptions);
        if (!is_wildcard) {
            auto query_type = QueryType::SPECIFIC_ELEMENT;
            auto atomic_number = get_value_for_equality_query(query);
            if (!is_atomic_number(atomic_number)) {
                throw UnrecognizedQueryError("Invalid atomic number");
            }
            auto element = Element(atomic_number);
            std::pair<QueryType, Element> new_vals = {query_type, element};
            std::pair<QueryType, Element> existing_vals = {
                query_props->query_type, query_props->element};
            check_value_for_conflicts(new_vals, existing_vals, QUERY_TYPE,
                                      seen_descriptions);
            query_props->query_type = query_type;
            query_props->element = element;
        }
    } else if (desc == "AtomFormalCharge") {
        throw_if_negated(query);
        auto val = get_value_for_equality_query(query);
        check_value_for_conflicts(val, query_props->charge, desc,
                                  seen_descriptions);
        query_props->charge = val;
    } else if (desc == "AtomHCount") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->total_h, desc,
                                  seen_descriptions);
        query_props->total_h = val;
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
        process_possible_wildcard_query(query, query_props, seen_descriptions);
        // TODO: this could probably also be an R-group - should check for
        //       R-group property before I call this function
    } else if (desc == "AtomRingBondCount") {
        throw_if_negated(query);
    } else if (desc == "AtomTotalDegree") {
        throw_if_negated(query);
        auto val = as_unsigned(get_value_for_equality_query(query));
        check_value_for_conflicts(val, query_props->num_connections, desc,
                                  seen_descriptions);
        query_props->num_connections = val;
    } else if (desc == "AtomType") {
        // AtomType combines both aromaticity and atomic number (aromaticity is
        // indicated by whether the value is over or under 1000)
        auto query_type = QueryType::SPECIFIC_ELEMENT;
        auto atom_type_num = get_value_for_equality_query(query);
        auto aromaticity = atom_type_num > ATOM_TYPE_AROMATIC_OFFSET
                               ? QueryAromaticity::AROMATIC
                               : QueryAromaticity::ALIPHATIC;
        auto atomic_number = atom_type_num % ATOM_TYPE_AROMATIC_OFFSET;
        if (!is_atomic_number(atomic_number)) {
            throw UnrecognizedQueryError("Invalid atomic number");
        }
        auto element = Element(atomic_number);
        std::pair<QueryType, Element> new_query_type_vals = {query_type,
                                                             element};
        std::pair<QueryType, Element> existing_query_type_vals = {
            query_props->query_type, query_props->element};
        check_value_for_conflicts(new_query_type_vals, existing_query_type_vals,
                                  QUERY_TYPE, seen_descriptions);
        check_value_for_conflicts(aromaticity, query_props->aromaticity,
                                  ATOM_AROMATICITY, seen_descriptions);
        query_props->query_type = query_type;
        query_props->element = element;
        query_props->aromaticity = aromaticity;
    } else {
        throw UnrecognizedQueryError("Description not recognized");
    }
}

/**
 * Read in all of the properties of the given query.
 */
static std::shared_ptr<AtomQueryProperties>
read_query(const RDKit::Atom::QUERYATOM_QUERY* const query)
{
    auto query_props = std::make_shared<AtomQueryProperties>();
    std::unordered_set<std::string> seen_descriptions;
    try {
        read_query_recursive(query, query_props, seen_descriptions);
    } catch (const UnrecognizedQueryError&) {
        // we couldn't parse the query, so throw out what we were able to
        // read (if anything) so that the failure is obvious
        query_props.reset(new AtomQueryProperties());
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
    if (!atom->hasQuery()) {
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
        // TODO: check for attachment point number
        props = read_query(query_atom->getQuery());
    }
    props->enhanced_stereo = read_enhanced_stereo_from_atom(atom);
    return props;
}

/**
 * Create a QueryAtom instance using the given properties. Also create any
 * queries necessary to describe the query type property.
 * @param query_props The query properties to use
 * @param queries A list of queries to be added to the returned QueryAtom. Any
 * queries needed to describe the query type property will be added to this
 * list, but *not* yet added to the atom.
 * @return A QueryAtom instance. This instance will be "empty" (i.e.
 * `QueryAtom()`) unless query_props specifies an R-group query type. In that
 * case, the returned atom will be an appropriately numbered attachment point
 * atom (since R-groups are specified differently from other queries).
 */
static std::shared_ptr<RDKit::Atom>
create_query_atom_and_query_for_query_type_property(
    const AtomQueryProperties* const query_props,
    std::vector<RDKit::Atom::QUERYATOM_QUERY*>& queries)
{
    std::shared_ptr<RDKit::Atom> query_atom;
    // first, create the atom
    if (query_props->query_type == QueryType::RGROUP) {
        query_atom = make_new_attachment_point(query_props->r_group);
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
            std::vector<RDKit::Atom::QUERYATOM_QUERY*> element_queries;
            for (auto element : query_props->allowed_list) {
                auto* cur_query =
                    RDKit::makeAtomNumQuery(static_cast<int>(element));
                element_queries.push_back(cur_query);
            }
            RDKit::Atom::QUERYATOM_QUERY* combined_query = nullptr;
            if (element_queries.size() == 1) {
                combined_query = *element_queries.begin();
            } else {
                combined_query = new RDKit::ATOM_OR_QUERY();
                for (auto* cur_query : element_queries) {
                    combined_query->addChild(
                        std::shared_ptr<RDKit::Atom::QUERYATOM_QUERY>(
                            cur_query));
                }
            }
            combined_query->setNegation(!allowed);
            queries.push_back(combined_query);
            break;
        }
        case QueryType::WILDCARD: {
            auto wildcard_query =
                ATOM_TOOL_QUERY_MAP.at(query_props->wildcard)();
            queries.push_back(wildcard_query);
            break;
        }
        case QueryType::SPECIFIC_ELEMENT: {
            auto atomic_num = static_cast<int>(query_props->element);
            query_atom->setAtomicNum(atomic_num);
            queries.push_back(RDKit::makeAtomNumQuery(atomic_num));
            break;
        }
        case QueryType::RGROUP: {
            // nothing to do here, since we already set the R-group when
            // creating the atom
            break;
        }
    }
    return query_atom;
}

/**
 * Update the atom and list of queries with all properties that were specified
 * in the General tab.
 * @param query_props The query properties to transfer to query_atom and queries
 * @param query_atom The query atom being created
 * @param queries A list of queries to be (eventually) added to query_atom. Any
 * queries created by this method will be added to this list, but *not* yet
 * added to the atom.
 */
static void create_queries_and_update_atom_for_general_properties(
    const AtomQueryProperties* const query_props,
    std::shared_ptr<RDKit::Atom> query_atom,
    std::vector<RDKit::Atom::QUERYATOM_QUERY*>& queries)
{
    if (query_props->isotope.has_value()) {
        queries.push_back(RDKit::makeAtomIsotopeQuery(*query_props->isotope));
        query_atom->setIsotope(*query_props->isotope);
    }
    if (query_props->charge.has_value()) {
        queries.push_back(
            RDKit::makeAtomFormalChargeQuery(*query_props->charge));
        query_atom->setFormalCharge(*query_props->charge);
    }
    if (query_props->unpaired_electrons.has_value()) {
        queries.push_back(RDKit::makeAtomNumRadicalElectronsQuery(
            *query_props->unpaired_electrons));
        query_atom->setNumRadicalElectrons(*query_props->unpaired_electrons);
    }
}

/**
 * Update the atom and list of queries with all properties that were specified
 * in the Advanced tab.
 * @param query_props The query properties to transfer to query_atom and queries
 * @param query_atom The query atom being created
 * @param queries A list of queries to be (eventually) added to query_atom. Any
 * queries created by this method will be added to this list, but *not* yet
 * added to the atom.
 */
static void create_queries_and_update_atom_for_advanced_properties(
    const AtomQueryProperties* const query_props,
    std::shared_ptr<RDKit::Atom> query_atom,
    std::vector<RDKit::Atom::QUERYATOM_QUERY*>& queries)
{
    if (query_props->total_h.has_value()) {
        queries.push_back(RDKit::makeAtomHCountQuery(*query_props->total_h));
        query_atom->setNumExplicitHs(*query_props->total_h);
    }
    if (query_props->num_connections.has_value()) {
        queries.push_back(
            RDKit::makeAtomTotalDegreeQuery(*query_props->num_connections));
    }

    if (query_props->aromaticity == QueryAromaticity::AROMATIC) {
        queries.push_back(RDKit::makeAtomAromaticQuery());
        query_atom->setIsAromatic(true);
    } else if (query_props->aromaticity == QueryAromaticity::ALIPHATIC) {
        queries.push_back(RDKit::makeAtomAliphaticQuery());
        query_atom->setIsAromatic(false);
    }

    if (query_props->ring_count_type == QueryCount::POSITIVE) {
        queries.push_back(RDKit::makeAtomInRingQuery());
    } else if (query_props->ring_count_type == QueryCount::EXACTLY) {
        queries.push_back(
            RDKit::makeAtomInNRingsQuery(query_props->ring_count_exact_val));
    }

    if (query_props->ring_bond_count_type == QueryCount::EXACTLY) {
        queries.push_back(RDKit::makeAtomRingBondCountQuery(
            query_props->ring_bond_count_exact_val));
    }
    if (query_props->smallest_ring_size.has_value()) {
        queries.push_back(
            RDKit::makeAtomMinRingSizeQuery(*query_props->smallest_ring_size));
    }

    auto smarts = boost::trim_copy(query_props->smarts_query);
    if (!smarts.empty()) {
        // TODO: need to store the SMARTS string somewhere on the atom, plus
        //       some way to skip trying to parse this query on read?
        auto smarts_atom = RDKit::SmartsToAtom(smarts);
        if (smarts_atom == nullptr) {
            for (auto cur_query : queries) {
                delete cur_query;
            }
            // we should never hit this since the dialog disables the OK
            // button when the SMARTS is invalid
            throw std::runtime_error("Invalid SMARTS pattern: " + smarts);
        }
        queries.push_back(smarts_atom->getQuery()->copy());
    }
}

// TODO: according to ChemicalKnowledge.cpp, product atoms can't have
//       allowed/disallowed list element queries, so it silently discards them.
//       (See SKETCH-843 and ENUM-402.)  Should probably handle this at the
//       export level.
/**
 * Create an RDKit atom with the given properties.
 */
std::shared_ptr<RDKit::Atom> create_atom_with_properties(
    const std::shared_ptr<AbstractAtomProperties> properties)
{
    // TODO: handle stereo - will need to return info about stereo group
    //       separately and pass it to MolModel, since it needs to be set on the
    //       molecule, not the atom
    if (!properties->isQuery()) {
        auto* atom_props = static_cast<const AtomProperties*>(properties.get());
        auto atom = std::make_shared<RDKit::Atom>(
            static_cast<int>(atom_props->element));
        atom->setIsotope(atom_props->isotope.value_or(0));
        atom->setFormalCharge(atom_props->charge);
        atom->setNumRadicalElectrons(atom_props->unpaired_electrons);
        return atom; //, atom_props.enhanced_stereo};
    } else {
        auto* query_props =
            static_cast<const AtomQueryProperties*>(properties.get());
        std::vector<RDKit::Atom::QUERYATOM_QUERY*> queries;
        auto query_atom = create_query_atom_and_query_for_query_type_property(
            query_props, queries);
        create_queries_and_update_atom_for_general_properties(
            query_props, query_atom, queries);
        create_queries_and_update_atom_for_advanced_properties(
            query_props, query_atom, queries);

        // add all of the queries to the atom
        for (auto* cur_query : queries) {
            if (!query_atom->hasQuery()) {
                query_atom->setQuery(cur_query);
            } else {
                query_atom->expandQuery(cur_query);
            }
        }
        return query_atom;
    }
}

} // namespace sketcher
} // namespace schrodinger
