#include "schrodinger/sketcher/rdkit/queries.h"

#include <array>
#include <memory>
#include <string_view>

#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/QueryBond.h>
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

namespace schrodinger
{
namespace sketcher
{

// Query type, negation allowed
constexpr std::array<std::pair<std::string_view, bool>, 12>
    supported_query_types{{
        {"AtomAtomicNum", true},
        {"AtomFormalCharge", false},
        {"AtomHCount", false},
        {"AtomInNRings", false},
        {"AtomInRing", true},
        {"AtomIsAliphatic", true},
        {"AtomIsAromatic", true},
        {"AtomMinRingSize", false},
        {"AtomNull", false},
        {"AtomRingBondCount", false},
        {"AtomTotalDegree", false},
        {"AtomType", true},
    }};

namespace
{
template <class T> bool get_query_values(
    const RDKit::Atom::QUERYATOM_QUERY& query, std::set<T>& value_collection,
    bool& is_negated,
    std::function<bool(const RDKit::Atom::QUERYATOM_QUERY&)> condition)
{
    const auto desc = query.getDescription();
    if (desc == "AtomAnd" || desc == "AtomOr") {
        for (auto it = query.beginChildren(); it != query.endChildren(); ++it) {
            if (get_query_values(**it, value_collection, is_negated,
                                 condition)) {
                is_negated |= query.getNegation();
            }
        }
    } else if (condition(query)) {
        // Not all queries we work on might be equality queries, but we just
        // want to get the value
        is_negated |= query.getNegation();
        auto equality_query =
            static_cast<const RDKit::ATOM_EQUALS_QUERY*>(&query);
        int value = equality_query->getVal();
        value_collection.insert(value);
        return true;
    }
    return false;
}

bool is_atom_type_query(const RDKit::Atom::QUERYATOM_QUERY& query)
{
    const auto desc = query.getDescription();
    if (desc == "AtomOr") {
        for (auto it = query.beginChildren(); it != query.endChildren(); ++it) {
            if (!is_atom_type_query(**it)) {
                return false;
            }
        }
    } else if (desc != "AtomAtomicNum" && desc != "AtomType") {
        return false;
    }
    return true;
}

} // namespace

bool is_query_fully_supported(const RDKit::Atom::QUERYATOM_QUERY& query)
{
    const auto desc = query.getDescription();
    if (desc == "AtomAnd") {
        for (auto it = query.beginChildren(); it != query.endChildren(); ++it) {
            if (!is_query_fully_supported(**it)) {
                return false;
            }
        }
    } else if (desc == "AtomOr") {
        // We only support OR queries for atomic numbers
        if (!is_atom_type_query(query)) {
            return false;
        }
    } else {
        auto check_query = [&query, &desc](const auto& query_type) {
            return desc == query_type.first &&
                   (query_type.second || !query.getNegation());
        };
        if (std::find_if(supported_query_types.begin(),
                         supported_query_types.end(),
                         check_query) == supported_query_types.end()) {
            return false;
        }
    }
    return true;
}

/* Insert into 'atomicNumbers' the atomic numbers of all chemical elements
 * included in the atom query. We ignore whether they are joined with
 * AND or OR, but in practice they should be OR since an actual atom can't
 * have more than one element! We also take negation on a global basis: if
 * one of the atomic numbers/types is negated, we assume all are, because we
 * cannot handle negations individually.
 */
void get_atomic_numbers_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                 std::set<unsigned>& atomic_numbers,
                                 bool& is_negated)
{
    std::set<unsigned> atom_types;
    auto is_atomic_num_query = [](const auto& query) {
        const auto desc = query.getDescription();
        return desc == "AtomAtomicNum" || desc == "AtomType";
    };
    get_query_values(query, atom_types, is_negated, is_atomic_num_query);

    for (auto atom_type : atom_types) {
        atomic_numbers.insert(atom_type % RDKIT_AROMATIC_ATOMIC_NUM_OFFSET);
    }
}

void get_formal_charges_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                 std::set<int>& charges)
{
    bool is_negated = false; // not used
    auto is_charge_query = [](const auto& query) {
        return !query.getNegation() &&
               query.getDescription() == "AtomFormalCharge";
    };
    get_query_values(query, charges, is_negated, is_charge_query);
}

void get_total_h_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                          std::set<unsigned>& total_hs)
{
    bool is_negated = false; // not used
    auto is_h_query = [](const auto& query) {
        return !query.getNegation() && query.getDescription() == "AtomHCount";
    };
    get_query_values(query, total_hs, is_negated, is_h_query);
}

void get_total_degree_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                               std::set<unsigned>& degree)
{
    bool is_negated = false; // not used
    auto is_degree_query = [](const auto& query) {
        return !query.getNegation() &&
               query.getDescription() == "AtomTotalDegree";
    };
    get_query_values(query, degree, is_negated, is_degree_query);
}

void get_num_rings_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                            std::set<unsigned>& num_rings)
{
    bool is_negated = false; // not used
    auto is_num_rings_query = [](const auto& query) {
        return !query.getNegation() && query.getDescription() == "AtomInNRings";
    };
    get_query_values(query, num_rings, is_negated, is_num_rings_query);
}

void get_in_ring_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                          std::set<unsigned>& in_ring, bool& is_negated)
{
    auto is_in_ring_query = [](const auto& query) {
        return query.getDescription() == "AtomInRing";
    };
    get_query_values(query, in_ring, is_negated, is_in_ring_query);
}

void get_ring_bond_count_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                  std::set<unsigned>& ring_bonds)
{
    bool is_negated = false; // not used
    auto is_ring_bond_count_query = [](const auto& query) {
        return !query.getNegation() &&
               query.getDescription() == "AtomRingBondCount";
    };
    get_query_values(query, ring_bonds, is_negated, is_ring_bond_count_query);
}

void get_min_ring_size_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                std::set<unsigned>& min_ring_size)
{
    bool is_negated = false; // not used
    auto is_min_ring_size_query = [](const auto& query) {
        return !query.getNegation() &&
               query.getDescription() == "AtomMinRingSize";
    };
    get_query_values(query, min_ring_size, is_negated, is_min_ring_size_query);
}

std::string get_atom_smarts(const RDKit::Atom* atom)
{
    return RDKit::SmartsWrite::GetAtomSmarts(atom);
}

void set_element_query(RDKit::Atom* atom, int atomic_number,
                       bool has_aromaticity, bool is_aromatic)
{
    atom->setAtomicNum(atomic_number);
    if (has_aromaticity) {
        atom->setIsAromatic(is_aromatic);
        atom->setQuery(RDKit::makeAtomTypeQuery(atomic_number, is_aromatic));
    } else {
        atom->setQuery(RDKit::makeAtomNumQuery(atomic_number));
    }
}

void set_element_list_query(RDKit::Atom* atom,
                            const std::set<unsigned>& atomic_numbers,
                            bool has_aromaticity, bool is_aromatic,
                            bool negated)
{
    for (auto atomicNumber : atomic_numbers) {
        RDKit::QueryAtom::QUERYATOM_QUERY* query;
        if (has_aromaticity) {
            query = RDKit::makeAtomTypeQuery(atomicNumber, is_aromatic);
        } else {
            query = RDKit::makeAtomNumQuery(atomicNumber);
        }

        if (atom->hasQuery()) {
            atom->expandQuery(query, Queries::COMPOSITE_OR);
        } else {
            atom->setAtomicNum(atomicNumber);
            atom->setQuery(query);
        }
    }
    atom->getQuery()->setNegation(negated);
}

void set_total_h_query(RDKit::Atom* atom, unsigned total_hs)
{
    auto query = RDKit::makeAtomHCountQuery(total_hs);
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_total_degree_query(RDKit::Atom* atom, unsigned degree)
{
    auto query = RDKit::makeAtomTotalDegreeQuery(degree);
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_num_rings_query(RDKit::Atom* atom, unsigned num_rings)
{
    auto query = RDKit::makeAtomInNRingsQuery(num_rings);
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_in_ring_query(RDKit::Atom* atom)
{
    auto query = RDKit::makeAtomInRingQuery();
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_ring_bond_count_query(RDKit::Atom* atom, unsigned ring_bonds)
{
    auto query = RDKit::makeAtomRingBondCountQuery(ring_bonds);
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_min_ring_size_query(RDKit::Atom* atom, unsigned min_ring_size)
{
    auto query = RDKit::makeAtomMinRingSizeQuery(min_ring_size);
    atom->expandQuery(query, Queries::COMPOSITE_AND);
}

void set_smarts_query(RDKit::Atom* atom, const std::string& smarts)
{
    std::unique_ptr<RDKit::Atom> smarts_atom{RDKit::SmartsToAtom(smarts)};
    if (smarts_atom != nullptr) {
        auto q_atom = static_cast<RDKit::QueryAtom*>(smarts_atom.get());
        atom->setQuery(q_atom->getQuery()->copy());
    }
}

} // namespace sketcher
} // namespace schrodinger
