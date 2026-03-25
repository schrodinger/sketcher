#include "schrodinger/rdkit_extensions/biologics_fingerprint.h"

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/MolOps.h>
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
    void update(std::string_view item)
    {
        size_t h1 = std::hash<std::string_view>{}(item);
        size_t h2 = boost::hash_value(item); // Or a different seed
        for (uint32_t i = 0; i < k_hashes; ++i) {
            // Kirsch-Mitzenmacher Optimization: h(i) = h1 + i * h2
            m_bitset.setBit((h1 + i * h2) % m_num_bits);
        }
    }
    void update(const fmt::memory_buffer& item)
    {
        update(std::string_view(item.data(), item.size()));
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

/**
 * @brief Appends unique atom key to buffer.
 *
 * Format: "monomer:{monomer_id}:{polymer_type}"
 * Does NOT clear the buffer - appends to existing content.
 */
void append_atom_key(const RDKit::Atom* atom, fmt::memory_buffer& buffer)
{
    std::string monomer_id;
    if (!atom->getPropIfPresent(ATOM_LABEL, monomer_id)) {
        throw std::invalid_argument(fmt::format("Atom {} missing '{}' property",
                                                atom->getIdx(), ATOM_LABEL));
    }

    const auto polymer_type = get_polymer_type(atom);
    fmt::format_to(std::back_inserter(buffer), "monomer:{}:{}", monomer_id,
                   polymer_type);
}

/**
 * @brief Appends unique bond key to buffer.
 *
 * Format: "bond:{atom1_key}:{atom2_key}:{bond_type}:{attachment}"
 * Does NOT clear the buffer - appends to existing content.
 */
void append_bond_key(const RDKit::Bond* bond, fmt::memory_buffer& buffer)
{
    std::string attachment_point;
    if (!bond->getPropIfPresent(LINKAGE, attachment_point)) {
        throw std::invalid_argument(fmt::format("Bond {} missing '{}' property",
                                                bond->getIdx(), LINKAGE));
    }

    const auto& mol = bond->getOwningMol();
    const auto* atom1 = mol.getAtomWithIdx(bond->getBeginAtomIdx());
    const auto* atom2 = mol.getAtomWithIdx(bond->getEndAtomIdx());

    fmt::format_to(std::back_inserter(buffer), "bond:");
    append_atom_key(atom1, buffer);
    fmt::format_to(std::back_inserter(buffer), ":");
    append_atom_key(atom2, buffer);
    fmt::format_to(std::back_inserter(buffer), ":{}:{}",
                   static_cast<int>(bond->getBondType()), attachment_point);
}

/// Finds the canonical rotation of a cyclic sequence using Booth's
/// algorithm.
int find_canonical_cycle_rotation(const std::vector<std::string>& atom_names)
{
    const auto n = static_cast<int>(atom_names.size());

    // Helper to access atom names with wraparound (treats sequence as circular)
    auto get_atom = [&](int idx) -> const std::string& {
        return atom_names[idx % n];
    };

    // Failure function for KMP-style matching
    std::vector<int> failure(n * 2, -1);

    // Candidate position for lexicographically minimal rotation
    int candidate_start = 0;

    // Compare positions in doubled sequence to find minimal rotation
    for (int pos = 1; pos < n * 2; ++pos) {
        int matched_len = failure[pos - candidate_start - 1];

        // Find longest prefix match using failure function
        while (matched_len != -1 &&
               get_atom(pos) != get_atom(candidate_start + matched_len + 1)) {
            if (get_atom(pos) < get_atom(candidate_start + matched_len + 1)) {
                candidate_start = pos - matched_len - 1;
            }
            matched_len = failure[matched_len];
        }

        // Update based on comparison result
        if (get_atom(pos) != get_atom(candidate_start + matched_len + 1)) {
            // matched_len must be -1 here (fell through from while loop)
            if (get_atom(pos) < get_atom(candidate_start + matched_len + 1)) {
                candidate_start = pos;
            }
            failure[pos - candidate_start] = -1;
        } else {
            // We have a match (atoms are equal)
            failure[pos - candidate_start] = matched_len + 1;
        }

        // Early termination once we've checked a full cycle
        if (pos - candidate_start >= n) {
            break;
        }
    }

    return candidate_start % n;
}

/**
 * @brief Appends unique cycle key to buffer.
 *
 * Format: "cycle:{size}:{bond1}:{bond2}:..."
 * Does NOT clear the buffer - appends to existing content.
 */
void append_cycle_key(const RDKit::ROMol& mol, const std::vector<int>& cycle,
                      fmt::memory_buffer& buffer)
{
    // Extract atom labels from the cycle
    std::vector<std::string> atom_names;
    atom_names.reserve(cycle.size());
    for (auto atom_idx : cycle) {
        auto atom = mol.getAtomWithIdx(atom_idx);
        atom_names.push_back(atom->getProp<std::string>(ATOM_LABEL));
    }

    // Build canonical key starting with cycle size
    fmt::format_to(std::back_inserter(buffer), "cycle:{}", cycle.size());

    // Find canonical starting position (lexicographically minimal rotation)
    auto start_idx = find_canonical_cycle_rotation(atom_names);

    // Append bond keys in canonical order
    for (auto i = 0u; i < cycle.size(); ++i) {
        fmt::format_to(std::back_inserter(buffer), ":");

        auto idx1 = cycle[(i + start_idx) % cycle.size()];
        auto idx2 = cycle[(i + start_idx + 1) % cycle.size()];
        auto bond = mol.getBondBetweenAtoms(idx1, idx2);
        append_bond_key(bond, buffer);
    }
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

void process_kmers(const RDKit::ROMol& mol, unsigned int k,
                   const std::function<void(const KmerPath&)>& callback)
{
    if (k < 2) {
        throw std::invalid_argument(
            fmt::format("k must be at least 2, got {}", k));
    }

    std::unordered_set<size_t> seen_paths; // Track unique paths
    KmerPath current_path;                 // Reused path vector
    current_path.reserve(k - 1);

    // Recursive DFS helper
    std::function<void(const RDKit::Atom*, const RDKit::Atom*)> dfs =
        [&](const RDKit::Atom* atom, const RDKit::Atom* parent) {
            // If we've reached target length, process the k-mer
            if (current_path.size() == k - 1) {
                const auto path_hash = hash_path(current_path);
                if (!seen_paths.contains(path_hash)) {
                    seen_paths.insert(path_hash);
                    callback(current_path);
                }
                return;
            }

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
                current_path.push_back(bond);
                dfs(neighbor, atom);
                current_path.pop_back(); // Backtrack
            }
        };

    // Start DFS from every monomer atom
    for (auto start_atom : mol.atoms()) {
        dfs(start_atom, nullptr);
    }
}

std::unique_ptr<ExplicitBitVect>
generate_biologics_fingerprint(const RDKit::ROMol& mol,
                               const BiologicsFingerprintConfig& config)
{
    // we want to limit the supported feature for now
    check_unsupported_features(mol);

    // this interface will be used to do the hashing
    BloomFilter bloom_filter{config.fp_size, config.num_hashes};

    // Reusable buffer for building keys
    fmt::memory_buffer key_buffer;

    // Add atom features
    for (auto atom : mol.atoms()) {
        key_buffer.clear();
        append_atom_key(atom, key_buffer);
        bloom_filter.update(key_buffer);
    }

    // Add k-mer features
    for (unsigned int k = 2; k <= config.max_k; ++k) {
        process_kmers(mol, k, [&](const KmerPath& kmer) {
            // Build composite key from all bonds in the k-mer
            key_buffer.clear();
            fmt::format_to(std::back_inserter(key_buffer), "sequence");
            for (auto bond : kmer) {
                fmt::format_to(std::back_inserter(key_buffer), ":");
                append_bond_key(bond, key_buffer);
            }
            bloom_filter.update(key_buffer);
        });
    }

    // Add cycle features
    constexpr bool includeDativeBonds = true;
    std::vector<std::vector<int>> rings;
    if (RDKit::MolOps::findSSSR(mol, rings, includeDativeBonds)) {
        for (const auto& cycle : rings) {
            key_buffer.clear();
            append_cycle_key(mol, cycle, key_buffer);
            bloom_filter.update(key_buffer);
        }
    }

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
