#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <rdkit/DataStructs/ExplicitBitVect.h>

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace RDKit
{
class ROMol;
class Atom;
} // namespace RDKit

namespace schrodinger::rdkit_extensions::fingerprint
{

/**
 * @brief Configuration for biologics fingerprint generation.
 */
struct BiologicsFingerprintConfig {
    size_t fp_size = 2048;       ///< Bit vector size
    unsigned int num_hashes = 3; ///< Hash functions per feature

    constexpr BiologicsFingerprintConfig() = default;
    constexpr BiologicsFingerprintConfig(size_t size, unsigned int hashes) :
        fp_size(size),
        num_hashes(hashes)
    {
    }
};

/**
 * @brief Generates a fingerprint for a biologics molecule.
 *
 * Creates a bit-vector fingerprint based on:
 * - Monomer composition (individual residues)
 *
 * NOTE: Fingerprints are count-less. Sequences differing only in repeat
 * count can have identical fingerprints when length exceeds max_k.
 * Example: "A×4" and "A×10" produce identical fingerprints when max_k=4.
 *
 * @param mol Biologics molecule with monomer annotations
 * @param config Biologics fingerprint configuration
 * @return Explicit bit vector fingerprint
 * @throws std::invalid_argument if molecule is not monomeric
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::unique_ptr<ExplicitBitVect>
generate_biologics_fingerprint(const RDKit::ROMol& mol,
                               const BiologicsFingerprintConfig& config = {});

/**
 * @brief Checks if query structure is a substructure of target.
 *
 * Converts HELM strings to molecules, generates fingerprints, and checks
 * if all bits set in the query fingerprint are also set in the target.
 *
 * @param query_helm Query structure in HELM notation
 * @param target_helm Target structure in HELM notation
 * @param config Biologics fingerprint configuration
 * @return True if query is a substructure of target
 * @throws std::invalid_argument if HELM parsing fails or molecules are not
 * monomeric
 */
[[nodiscard]] RDKIT_EXTENSIONS_API bool is_substructure_fingerprint_match(
    const std::string& query_helm, const std::string& target_helm,
    const BiologicsFingerprintConfig& config = {});

} // namespace schrodinger::rdkit_extensions::fingerprint
