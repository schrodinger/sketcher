#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/token_scanner.h"
#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"

namespace helm
{

struct polymer_group {
    std::string_view id;
    std::string_view items;
    bool is_polymer_union;
};

struct connection {
    std::string_view from_id;
    std::string_view from_res;
    std::string_view from_rgroup;
    std::string_view to_id;
    std::string_view to_res;
    std::string_view to_rgroup;
    std::string_view annotation;
};

struct monomer {
    std::string_view id;
    bool is_smiles = false;
    bool is_branch = false;
    bool is_list = false;
    std::string_view annotation;
};

struct repetition {
    size_t start;
    size_t size;
    std::string_view num_repetitions;
    std::string_view annotation;
};

struct polymer {
    std::string_view id;
    std::vector<monomer> monomers;
    std::vector<repetition> repetitions;
    std::unordered_set<size_t> wildcard_and_unknown_residues;
    std::unordered_set<std::string_view> residue_names;
    std::string_view annotation;
};

struct helm_info {
    std::vector<polymer> polymers;
    std::vector<connection> connections;
    std::vector<polymer_group> polymer_groups;
    std::string_view extended_annotations;
};

class RDKIT_EXTENSIONS_API HelmParser
{
  public:
    HelmParser(const std::string_view input_helm) : m_input(input_helm){};
    // top-level function to initiate parsing
    helm_info parse();

    // apis for polymers
    void add_monomer(const std::string_view monomer_id, const bool is_smiles,
                     const bool is_branch, const bool is_list,
                     const std::string_view annotation);

    void add_monomer_with_id(const std::string_view monomer_id,
                             const std::string_view annotation);
    void add_smiles_monomer(const std::string_view monomer_id,
                            const std::string_view annotation);
    void add_monomer_list(const std::string_view monomer_id,
                          const std::string_view annotation);
    void add_residue_name(const std::string_view monomer_id);
    void add_wildcard_or_unknown_residue();

    void mark_branch_monomer(const size_t branch_group_size);
    void
    mark_last_n_monomers_as_repeated(const size_t repetition_size,
                                     const std::string_view num_repetitions,
                                     const std::string_view annotation);
    void add_polymer(const std::string_view polymer_id,
                     const std::string_view annotation);

    // apis for connections
    void add_connection(connection&& connection);

    // apis for polymer groups
    void add_polymer_group(const std::string_view polymer_group_id,
                           const std::string_view items,
                           const bool is_polymer_union);

    // apis for extended annotations
    void add_extended_annotations(const std::string_view extended_annotations);

    bool hasErrors();
    void saveError(const std::string_view& failed_token,
                   const std::string& err_msg);
    void saveError(const unsigned int num_chars_processed,
                   const std::string& err_msg);

  private:
    std::string_view m_input;
    std::vector<std::string> m_errors;
    helm_info m_parsed_info;
};

} /* end namespace helm */
