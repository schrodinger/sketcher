#include "schrodinger/rdkit_extensions/biologics_fingerprint.h"

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/AtomIterators.h>

#include <boost/container_hash/hash.hpp>
#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

namespace schrodinger::rdkit_extensions::fingerprint
{

namespace
{

/**
 * @brief Bloom Filter interface for hashing fingerprint features.
 * * PERFORMANCE NOTES:
 * - Uses Hybrid Hashing (1 heavy hash + k light salts).
 */
class BloomFilter
{
  private:
    ExplicitBitVect m_bitset;
    const size_t m_num_bits;
    uint32_t k_hashes;

    // 64-bit Golden Ratio constant for bit distribution.
    static constexpr uint64_t GOLDEN_RATIO_64 = 0x9e3779b97f4a7c15ULL;

  public:
    /**
     * @brief Constructor
     * @param m Total bits (suggested: 10 * expected_elements for 1% FPR).
     * @param k Number of hashes (suggested: 7 for 1% FPR).
     */
    explicit BloomFilter(size_t m, uint32_t k) :
        m_bitset(m),
        m_num_bits(m),
        k_hashes(k)
    {
    }

    /**
     * @brief High-throughput insertion.
     * Hashes the object once, then perturbs the seed k times.
     */
    void update(const std::string& item)
    {
        const size_t primary_hash = std::hash<std::string>{}(item);

        for (uint32_t i = 0; i < k_hashes; ++i) {
            size_t seed = primary_hash;

            // Salted perturbation using the iteration index and Golden Ratio.
            size_t salt = i * GOLDEN_RATIO_64;
            boost::hash_combine(seed, salt);

            m_bitset.setBit(seed % m_num_bits);
        }
    }

    /**
     * @brief API to retrieve the underlying bitset.
     * Allows for direct serialization or bitwise operations between filters.
     */
    std::unique_ptr<ExplicitBitVect> get_bitset() const noexcept
    {
        return std::make_unique<ExplicitBitVect>(m_bitset);
    }
};

/// Gets biologics type ID without trailing digits
std::string get_polymer_type(const RDKit::Atom* atom)
{
    // Get biologics type ID and strip trailing digits
    std::string polymer_id = get_polymer_id(atom);
    return polymer_id.substr(0, polymer_id.find_first_of("0123456789"));
}

/// Generates unique key for an atom
std::string get_atom_key(const RDKit::Atom* atom)
{
    std::string monomer_id;
    if (!atom->getPropIfPresent(ATOM_LABEL, monomer_id)) {
        throw std::invalid_argument(fmt::format("Atom {} missing '{}' property",
                                                atom->getIdx(), ATOM_LABEL));
    }

    const auto polymer_type = get_polymer_type(atom);
    return fmt::format("monomer:{}:{}", monomer_id, polymer_type);
}

void check_unsupported_features(const RDKit::ROMol& mol)
{
    // Verify molecule is monomeric
    auto atoms_it = mol.atoms();
    std::vector<RDKit::Atom*> atoms{atoms_it.begin(), atoms_it.end()};
    if (!isMonomeric(mol) || !std::ranges::any_of(atoms, [](auto atom) {
            auto info = atom->getMonomerInfo();
            return info && static_cast<const RDKit::AtomPDBResidueInfo*>(info);
        })) {
        throw std::invalid_argument("Molecule must be monomeric (biologics)");
    }

    // monomer lists are not supported
    if (std::ranges::any_of(
            atoms, [](auto atom) { return atom->hasProp(MONOMER_LIST); })) {
        throw std::invalid_argument(
            "Biologics fingerprinting doesn't support monomer lists.");
    }

    // repetitions are not supported
    if (std::ranges::any_of(atoms, [](auto atom) {
            return atom->hasProp(REPETITION_DUMMY_ID);
        })) {
        throw std::invalid_argument(
            "Biologics fingerprinting doesn't support monomer repetitions.");
    }

    // unknown monomers are not supported
    if (std::ranges::any_of(
            atoms, [](auto atom) { return atom->hasProp(UNKNOWN_MONOMER); })) {
        throw std::invalid_argument(
            "Biologics fingerprinting doesn't support unknown monomers.");
    }
}

} // anonymous namespace

std::unique_ptr<ExplicitBitVect>
generate_biologics_fingerprint(const RDKit::ROMol& mol,
                               const BiologicsFingerprintConfig& config)
{
    // we want to limit the supported feature for now
    check_unsupported_features(mol);

    // this interface will be used to do the hashing
    BloomFilter bloom_filter{config.fp_size, config.num_hashes};

    // Add atom features
    for (auto atom : mol.atoms()) {
        const auto atom_key = get_atom_key(atom);
        bloom_filter.update(atom_key);
    }

    // TODO: Add k-mer features

    // TODO: Add cycle features

    return bloom_filter.get_bitset();
}

bool is_substructure_fingerprint_match(const std::string& query_helm,
                                       const std::string& target_helm,
                                       const BiologicsFingerprintConfig& config)
{
    // Convert HELM to RDKit molecules
    auto query_mol = helm::helm_to_rdkit(query_helm);
    if (!query_mol) {
        throw std::invalid_argument(
            fmt::format("Failed to parse query HELM: {}", query_helm));
    }

    auto target_mol = helm::helm_to_rdkit(target_helm);
    if (!target_mol) {
        throw std::invalid_argument(
            fmt::format("Failed to parse target HELM: {}", target_helm));
    }

    // Generate fingerprints
    auto query_fp = generate_biologics_fingerprint(*query_mol, config);
    auto target_fp = generate_biologics_fingerprint(*target_mol, config);

    // Check subset relationship
    return (*query_fp & *target_fp) == *query_fp;
}

} // namespace schrodinger::rdkit_extensions::fingerprint
