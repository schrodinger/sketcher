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
        // if the label contains only one of the above characters, either alone
        // or as an AND query with something else, we can use it as a bond type
        QString query_label_qstring = QString::fromStdString(query_label);
        auto list = query_label_qstring.split("+");
        auto result_itr = query_label_to_bond_type.end();
        for (auto item : list) {
            auto bond_type_itr =
                query_label_to_bond_type.find(item.toStdString());
            if (bond_type_itr != query_label_to_bond_type.end()) {
                if (result_itr != query_label_to_bond_type.end()) {
                    // more than one bond type found
                    result_itr = query_label_to_bond_type.end();
                    break;
                } else {
                    result_itr = bond_type_itr;
                }
            }
        }
        if (result_itr != query_label_to_bond_type.end()) {
            bond_type = result_itr->second;
            list.removeAll(result_itr->first.c_str());
            query_label = list.join("+").toStdString();
        } else {
            bond_type = RDKit::Bond::BondType::SINGLE;
        }
    }
    return {bond_type, query_label};
}

BondTopology get_bond_topology(const RDKit::Bond* const bond)
{
    if (!bond->hasQuery()) {
        return BondTopology::EITHER;
    }
    auto query = bond->getQuery();
    for (auto child = query->beginChildren(); child != query->endChildren();
         ++child) {
        if ((*child)->getDescription() == "BondInRing") {
            return ((*child)->getNegation() ? BondTopology::NOT_IN_RING
                                            : BondTopology::IN_RING);
        }
    }
    return BondTopology::EITHER;
}

void set_bond_topology(RDKit::QueryBond* const bond, BondTopology topology)
{
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> query(
        RDKit::makeBondIsInRingQuery());
    query->setNegation(topology == BondTopology::NOT_IN_RING);
    if (bond->hasQuery()) {
        // if topology is already set, replace it
        for (auto child = bond->getQuery()->beginChildren();
             child != bond->getQuery()->endChildren(); ++child) {
            if ((*child)->getDescription() == query->getDescription()) {
                (*child)->setNegation(query->getNegation());
                return;
            }
        }
        // otherwise add it
        bond->expandQuery(query->copy(), Queries::COMPOSITE_AND);
    } else {
        bond->setQuery(query->copy());
    }
}

std::shared_ptr<RDKit::Bond>
make_new_bond_without_topology(const RDKit::QueryBond* existing_bond)
{
    auto bond = std::make_shared<RDKit::QueryBond>(*existing_bond);
    // we are assuming that the existing bond is an AND query between two
    // children.
    if (!existing_bond->hasQuery()) {
        return bond;
    }
    auto existing_query = existing_bond->getQuery();
    if (existing_query->getDescription() != "BondAnd") {

        return bond;
    }

    int child_count = 0;
    auto child_to_keep_itr = existing_query->endChildren();
    for (auto ci = existing_query->beginChildren();
         ci != existing_query->endChildren(); ++ci) {
        ++child_count;
        if ((*ci)->getDescription() != "BondInRing") {
            // we want to keep this child, but only if there's exactly one
            //  non-BondInRing child and a total of two children (i.e. the other
            //  child is a BondInRing query)
            if (child_to_keep_itr == existing_query->endChildren()) {
                child_to_keep_itr = ci;
            } else {
                child_to_keep_itr = existing_query->endChildren();
            }
        }
    }
    if (child_count != 2 ||
        child_to_keep_itr == existing_query->endChildren()) {
        return bond;
    }
    // if the only remaining child is a BondOrder query, we return a normal bond
    // of the corresponding order
    if ((*child_to_keep_itr)->getDescription() == "BondOrder") {
        auto equals_query = dynamic_cast<const RDKit::BOND_EQUALS_QUERY*>(
            child_to_keep_itr->get());
        if (equals_query == nullptr) {
            return bond;
        }
        switch (equals_query->getVal()) {
            case 1:
                return std::make_shared<RDKit::Bond>(
                    RDKit::Bond::BondType::SINGLE);
            case 2:
                return std::make_shared<RDKit::Bond>(
                    RDKit::Bond::BondType::DOUBLE);
            case 3:
                return std::make_shared<RDKit::Bond>(
                    RDKit::Bond::BondType::TRIPLE);
            case 12:
                return std::make_shared<RDKit::Bond>(
                    RDKit::Bond::BondType::AROMATIC);
            default:
                return bond;
        }
    }
    // otherwise we return a query bond with the remaining child
    if (child_count == 2 &&
        child_to_keep_itr != existing_query->endChildren()) {
        bond->setQuery(child_to_keep_itr->get()->copy());
    }
    return bond;
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
    auto desc_label = parse_description_for_bond_query_node(query);

    if (query->getNegation()) {
        // Use a specific symbol for BondTopology::NOT_IN_RING
        if (query->getDescription() == "BondInRing") {
            return "Not ⭔";
        }
        // Otherwise use bang for negation
        return "!" + desc_label;
    }

    return desc_label;
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
        return "⭔";
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
