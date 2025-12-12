#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

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
    std::string_view annotation = {};
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

/**
 * @brief Constructs a detailed, user-friendly error message, including a
 *        truncated snippet of the input string and an error pointer.
 *
 * @param err_msg The high-level error description
 * @param pos The 0-based index where the error occurred.
 * @param input The fill input string
 * @return A formatted error string suitable for output
 */
[[nodiscard]] std::string construct_error_msg(const std::string_view err_msg,
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
construct_error_msg(const std::string_view err_msg,
                    const std::string_view& failed_token,
                    const std::string_view& input);

/// Helper api to determine whether a multi-character monomer token is a valid
/// inline SMILES monomer.
[[nodiscard]] bool is_smiles_monomer(const std::string_view&);

/**
 * @brief Parses a HELM v2 string into a structured representation.
 *
 * Produces a `helm_info` object describing polymers, connections,
 * annotations, and version fields.
 *
 * @param input  HELM v2 input string (e.g. "BLOB1{Bead}$$$$V2.0").
 * @param do_throw  Throw on parse errors if true; otherwise log and return
 *        std::nullopt (default: true).
 *
 * @return Parsed HELM info, or std::nullopt on failure.
 * @throws std::invalid_argument if parsing fails and `doThrow == true`.
 */
RDKIT_EXTENSIONS_API std::optional<helm_info> parse_helm(std::string_view input,
                                                         bool do_throw = true);
} /* end namespace helm */
