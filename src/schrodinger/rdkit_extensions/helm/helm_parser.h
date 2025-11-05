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
    std::string_view helm_version;
};

class RDKIT_EXTENSIONS_API HelmParser
{
  public:
    HelmParser(const std::string_view input_helm);
    // top-level function to initiate parsing
    helm_info parse();

    // apis for polymers
    void addMonomer(const std::string_view monomer_id, const bool is_smiles,
                    const bool is_branch, const bool is_list,
                    const std::string_view annotation);

    void addMonomerWithId(const std::string_view monomer_id,
                          const std::string_view annotation);
    void addSmilesMonomer(const std::string_view monomer_id,
                          const std::string_view annotation);
    void addMonomerList(const std::string_view monomer_id,
                        const std::string_view annotation);
    void addResidueName(const std::string_view monomer_id);
    void addWildcardOrUnknownResidue();

    void markBranchMonomer(const size_t branch_group_size);
    void markLastNMonomersAsRepeated(const size_t repetition_size,
                                     const std::string_view num_repetitions,
                                     const std::string_view annotation);
    void addPolymer(const std::string_view polymer_id,
                    const std::string_view annotation);

    // apis for connections
    void addConnection(connection&& connection);

    // apis for polymer groups
    void addPolymerGroup(const std::string_view polymer_group_id,
                         const std::string_view items,
                         const bool is_polymer_union);

    // apis for extended annotations
    void addExtendedAnnotations(const std::string_view extended_annotations);

    // apis for version
    void addHelmVersion(const std::string_view helm_version);

    void saveError(const std::string_view& failed_token,
                   const std::string& err_msg);
    void saveError(const unsigned int num_chars_processed,
                   const std::string& err_msg);

  private:
    std::string_view m_input;
    std::vector<std::string> m_errors;
    helm_info m_parsed_info;
};

/**
 * @brief Constructs a detailed, user-friendly error message, including a
 *        truncated snippet of the input string and an error pointer.
 *
 * @param err_msg The high-level error description
 * @param pos The 0-based index where the error occurred.
 * @param input The fill input string
 * @return A formatted error string suitable for output
 */
[[nodiscard]] std::string construct_error_msg(const std::string& err_msg,
                                              const unsigned int pos,
                                              const std::string_view& input);
/**
 * @brief Constructs a detailed, user-friendly error message, including a
 *        truncated snippet of the input string and an error pointer.
 *
 * @param err_msg The high-level error description
 * @param failed_token The rejected token. This must point to a location within
 *        the input value i.e. the memory address of the underlying value must
 *        fall under that of the input string.
 * @param input The fill input string
 * @return A formatted error string suitable for output
 */
[[nodiscard]] std::string
construct_error_msg(const std::string& err_msg,
                    const std::string_view& failed_token,
                    const std::string_view& input);

/// Helper api to determine whether a multi-character monomer token is a valid
/// inline SMILES monomer.
[[nodiscard]] bool is_smiles_monomer(const std::string_view&);

} /* end namespace helm */
