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

/// DFS state: (current_atom, parent_atom, bond_path)
struct DFSState {
    const RDKit::Atom* atom;
    const RDKit::Atom* parent;
    KmerPath path;
};

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
        size_t h1 = std::hash<std::string>{}(item);
        size_t h2 = boost::hash_value(item); // Or a different seed
        for (uint32_t i = 0; i < k_hashes; ++i) {
            // Kirsch-Mitzenmacher Optimization: h(i) = h1 + i * h2
            m_bitset.setBit((h1 + i * h2) % m_num_bits);
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

/// Generates unique key for a bond
std::string get_bond_key(const RDKit::Bond* bond)
{
    std::string attachment_point;
    if (!bond->getPropIfPresent(LINKAGE, attachment_point)) {
        throw std::invalid_argument(fmt::format("Bond {} missing '{}' property",
                                                bond->getIdx(), LINKAGE));
    }

    const auto& mol = bond->getOwningMol();
    const auto* atom1 = mol.getAtomWithIdx(bond->getBeginAtomIdx());
    const auto* atom2 = mol.getAtomWithIdx(bond->getEndAtomIdx());

    const auto bond_type = bond->getBondType();
    return fmt::format("bond:{}:{}:{}:{}", get_atom_key(atom1),
                       get_atom_key(atom2), static_cast<int>(bond_type),
                       attachment_point);
}

/// Creates unique identifier for a bond path
size_t hash_path(const KmerPath& path)
{
    size_t seed = 0;
    for (const auto* bond : path) {
        // Combine bond indices using simple hash combine
        boost::hash_combine(seed, bond->getIdx());
    }
    return seed;
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

/// Checks if bond is a valid biologics bond
bool is_interesting_bond(const RDKit::Bond* bond)
{
    std::string attachment_point;
    // Exclude hydrogen bonds
    return bond->getPropIfPresent(LINKAGE, attachment_point) &&
           !attachment_point.starts_with("pair");
}

} // anonymous namespace
  //
std::vector<KmerPath> extract_kmers(const RDKit::ROMol& mol, unsigned int k)
{
    if (k < 2) {
        throw std::invalid_argument(
            fmt::format("k must be at least 2, got {}", k));
    }

    std::vector<KmerPath> paths;
    std::unordered_set<size_t> seen_paths; // Track unique paths

    // Start DFS from every monomer atom
    for (auto start_atom : mol.atoms()) {
        std::vector<DFSState> stack;
        stack.push_back({start_atom, nullptr, {}});

        while (!stack.empty()) {
            auto [atom, parent, path] = std::move(stack.back());
            stack.pop_back();

            // Explore neighbors via biologics bonds
            for (const auto* bond : mol.atomBonds(atom)) {
                if (!is_interesting_bond(bond)) {
                    continue;
                }

                // Only traverse bonds in the forward direction
                if (bond->getBeginAtomIdx() != atom->getIdx()) {
                    continue;
                }

                const auto* neighbor =
                    mol.getAtomWithIdx(bond->getEndAtomIdx());

                // Don't go backwards
                if (neighbor == parent) {
                    continue;
                }

                // Extend the path
                KmerPath new_path = path;
                new_path.push_back(bond);

                // If we've reached the target k-mer length, record it
                if (new_path.size() == k - 1) {
                    const auto path_hash = hash_path(new_path);
                    if (!seen_paths.contains(path_hash)) {
                        seen_paths.insert(path_hash);
                        paths.push_back(std::move(new_path));
                    }
                } else {
                    stack.push_back({neighbor, atom, std::move(new_path)});
                }
            }
        }
    }

    return paths;
}

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

    // Add k-mer features
    for (unsigned int k = 2; k <= config.max_k; ++k) {
        std::ranges::for_each(extract_kmers(mol, k), [&](auto& kmer) {
            // Build composite key from all bonds in the k-mer
            std::string kmer_key = "sequence";
            for (auto bond : kmer) {
                kmer_key += ":";
                kmer_key += get_bond_key(bond);
            }

            bloom_filter.update(kmer_key);
        });
    }

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
