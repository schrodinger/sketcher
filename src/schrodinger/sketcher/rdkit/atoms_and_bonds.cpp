#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

#include <algorithm>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

namespace schrodinger
{
namespace sketcher
{

namespace
{

class UnsupportedQueryException : public std::exception
{
};

} // namespace

static std::string
parse_bond_query_node(const RDKit::Bond::QUERYBOND_QUERY* const query);

static std::string parse_all_children_of_bond_query_node(
    const RDKit::Bond::QUERYBOND_QUERY* const query, std::string join_string);

static std::string parse_description_for_bond_query_node(
    const RDKit::Bond::QUERYBOND_QUERY* const query);

std::pair<RDKit::Bond::BondType, std::string>
get_bond_type_and_query_label(const RDKit::Bond* const bond)
{
    RDKit::Bond::BondType bond_type;
    std::string query_label = "";
    if (!bond->hasQuery()) {
        bond_type = bond->getBondType();
    } else {
        query_label = get_label_for_bond_query(bond->getQuery());
        static const std::unordered_map<std::string, RDKit::Bond::BondType>
            query_label_to_bond_type{
                {"S", RDKit::Bond::BondType::SINGLE},
                {"D", RDKit::Bond::BondType::DOUBLE},
                {"T", RDKit::Bond::BondType::TRIPLE},
                {"A", RDKit::Bond::BondType::AROMATIC},
            };
        auto bond_type_itr = query_label_to_bond_type.find(query_label);
        if (bond_type_itr != query_label_to_bond_type.end()) {
            // if the query specifies exactly one bond type and nothing more,
            // then render the bond as that type instead of displaying a query
            // label
            bond_type = bond_type_itr->second;
            query_label = "";
        } else {
            // If a bond has a query label, always draw it as a single bond
            bond_type = RDKit::Bond::BondType::SINGLE;
        }
    }
    return {bond_type, query_label};
}

std::string
get_label_for_bond_query(const RDKit::Bond::QUERYBOND_QUERY* const query)
{
    try {
        return parse_bond_query_node(query);
    } catch (const UnsupportedQueryException&) {
        // Display "Query" for all queries that we can't parse
        return "Query";
    }
}

/**
 * Get the label string to display for the given bond query
 *
 * @param query The query to parse.  Must be non-null.
 * @throw UnsupportedQueryException if the query cannot be parsed
 */
static std::string
parse_bond_query_node(const RDKit::Bond::QUERYBOND_QUERY* const query)
{
    std::string negation;
    if (query->getNegation()) {
        negation = "!";
    }
    return negation + parse_description_for_bond_query_node(query);
}

/**
 * Get the label string to display for the description of the given bond query.
 * Note that this function does *not* account for any negation of the query.
 *
 * @param query The query to parse.  Must be non-null.
 * @throw UnsupportedQueryException if the query cannot be parsed
 */
static std::string parse_description_for_bond_query_node(
    const RDKit::Bond::QUERYBOND_QUERY* const query)
{
    auto desc = query->getDescription();
    if (desc == "BondNull") {
        return "Any";
    } else if (desc == "SingleOrDoubleBond") {
        return "S/D";
    } else if (desc == "SingleOrAromaticBond") {
        return "S/A";
    } else if (desc == "DoubleOrAromaticBond") {
        return "D/A";
    } else if (desc == "SingleOrDoubleOrAromaticBond") {
        return "S/D/A";
    } else if (desc == "BondInRing") {
        return "Ring";
    } else if (desc == "BondOrder") {
        auto equals_query =
            dynamic_cast<const RDKit::BOND_EQUALS_QUERY*>(query);
        if (equals_query == nullptr) {
            throw UnsupportedQueryException();
        }
        switch (equals_query->getVal()) {
            // these values were determined by converting SMARTS queries to
            // RDKit objects and introspecting the resulting bond queries
            case 1:
                // single bond
                return "S";
            case 2:
                // double bond
                return "D";
            case 3:
                // triple bond
                return "T";
            case 12:
                // aromatic bond - the twelve is intended as a one and a two
                // (i.e. a list of two single digit numbers) since an aromatic
                // bond is "between" a single bond and a double bond
                return "A";
            default:
                throw UnsupportedQueryException();
        }
    } else if (desc == "BondOr") {
        return parse_all_children_of_bond_query_node(query, "/");
    } else if (desc == "BondAnd") {
        return parse_all_children_of_bond_query_node(query, "+");
    }
    throw UnsupportedQueryException();
}

/**
 * Generate label strings for all children of the given bond query and join
 * those descriptions together using join_string.
 *
 * @param query The query to parse.  Must be non-null.
 * @param join_string The string to put in between each child description
 * @throw UnsupportedQueryException if the query cannot be parsed
 */
static std::string parse_all_children_of_bond_query_node(
    const RDKit::Bond::QUERYBOND_QUERY* const query, std::string join_string)
{
    std::vector<std::string> child_text;
    std::transform(query->beginChildren(), query->endChildren(),
                   std::back_inserter(child_text), [](const auto ci) {
                       return parse_bond_query_node(ci.get());
                   });
    return boost::algorithm::join(child_text, join_string);
}

bool has_any_implicit_Hs(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    // RDKit calls some hydrogens "explicit" but they are not really
    // explicit in the sense of being in the graph and we need to count them
    // in. We count all hydrogens that are not in the graph as implicit.
    return std::any_of(atoms.begin(), atoms.end(),
                       [](auto atom) { return atom->getTotalNumHs() > 0; });
}

} // namespace sketcher
} // namespace schrodinger
