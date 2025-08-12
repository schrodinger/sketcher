
#pragma once

#include <set>
#include <string>

#include <rdkit/GraphMol/Atom.h>

namespace schrodinger
{
namespace sketcher
{

constexpr int RDKIT_AROMATIC_ATOMIC_NUM_OFFSET = 1000;

bool is_query_fully_supported(const RDKit::Atom::QUERYATOM_QUERY& query);

void get_atomic_numbers_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                 std::set<unsigned>& atomic_numbers,
                                 bool& is_negated);

void get_formal_charges_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                 std::set<int>& charges);

void get_total_h_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                          std::set<unsigned>& total_hs);

void get_total_degree_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                               std::set<unsigned>& degree);

void get_num_rings_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                            std::set<unsigned>& num_rings);

void get_in_ring_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                          std::set<unsigned>& in_ring, bool& is_negated);

void get_ring_bond_count_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                  std::set<unsigned>& ring_bonds);

void get_min_ring_size_of_query(const RDKit::Atom::QUERYATOM_QUERY& query,
                                std::set<unsigned>& min_ring_size);

std::string get_atom_smarts(const RDKit::Atom* atom);

void set_element_query(RDKit::Atom* atom, int atomic_number,
                       bool has_aromaticity, bool is_aromatic);

void set_element_list_query(RDKit::Atom* atom,
                            const std::set<unsigned>& atomic_numbers,
                            bool has_aromaticity, bool is_aromatic,
                            bool negated);

void set_total_h_query(RDKit::Atom* atom, unsigned total_hs);

void set_total_degree_query(RDKit::Atom* atom, unsigned degree);

void set_num_rings_query(RDKit::Atom* atom, unsigned num_rings);

void set_in_ring_query(RDKit::Atom* atom);

void set_ring_bond_count_query(RDKit::Atom* atom, unsigned ring_bonds);

void set_min_ring_size_query(RDKit::Atom* atom, unsigned min_ring_size);

void set_smarts_query(RDKit::Atom* atom, const std::string& smarts);

} // namespace sketcher
} // namespace schrodinger
