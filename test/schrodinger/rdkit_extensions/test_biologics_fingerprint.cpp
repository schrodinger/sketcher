#define BOOST_TEST_MODULE test_biologics_fingerprint

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/rdkit_extensions/biologics_fingerprint.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

using namespace schrodinger::rdkit_extensions::fingerprint;
namespace bdata = boost::unit_test::data;

namespace
{

/// Helper to generate fingerprint from HELM string
auto generate_fp(const std::string& helm_string,
                 const BiologicsFingerprintConfig& config = {})
{
    auto mol = helm::helm_to_rdkit(helm_string);
    return generate_biologics_fingerprint(*mol, config);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(kmer_extraction)

BOOST_AUTO_TEST_CASE(linear_peptide_kmer_counts)
{
    // For linear "A.C.G.T":
    auto mol = helm::helm_to_rdkit("PEPTIDE1{A.C.G.T}$$$$V2.0");

    // k=2: A-C, C-G, G-T (3 k-mers)
    std::vector<KmerPath> kmers_k2;
    process_kmers(*mol, 2,
                  [&](const KmerPath& kmer) { kmers_k2.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k2.size(), 3);

    // k=3: A-C-G, C-G-T (2 k-mers)
    std::vector<KmerPath> kmers_k3;
    process_kmers(*mol, 3,
                  [&](const KmerPath& kmer) { kmers_k3.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k3.size(), 2);

    // k=4: A-C-G-T (1 k-mer)
    std::vector<KmerPath> kmers_k4;
    process_kmers(*mol, 4,
                  [&](const KmerPath& kmer) { kmers_k4.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k4.size(), 1);

    // k=5: (0 k-mers)
    std::vector<KmerPath> kmers_k5;
    process_kmers(*mol, 5,
                  [&](const KmerPath& kmer) { kmers_k5.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k5.size(), 0);
}

BOOST_AUTO_TEST_CASE(longer_sequence_kmer_counts)
{
    // For 6-mer: A.C.G.T.E.F
    auto mol = helm::helm_to_rdkit("PEPTIDE1{A.C.G.T.E.F}$$$$V2.0");

    // k=2: 5 k-mers
    std::vector<KmerPath> kmers_k2;
    process_kmers(*mol, 2,
                  [&](const KmerPath& kmer) { kmers_k2.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k2.size(), 5);

    // k=3: 4 k-mers
    std::vector<KmerPath> kmers_k3;
    process_kmers(*mol, 3,
                  [&](const KmerPath& kmer) { kmers_k3.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k3.size(), 4);

    // k=4: 3 k-mers
    std::vector<KmerPath> kmers_k4;
    process_kmers(*mol, 4,
                  [&](const KmerPath& kmer) { kmers_k4.push_back(kmer); });
    BOOST_CHECK_EQUAL(kmers_k4.size(), 3);
}

BOOST_AUTO_TEST_CASE(kmer_path_correctness)
{
    auto mol = helm::helm_to_rdkit("PEPTIDE1{A.C.G}$$$$V2.0");

    std::vector<KmerPath> kmers;
    process_kmers(*mol, 2,
                  [&](const KmerPath& kmer) { kmers.push_back(kmer); });

    // Each k=2 path should have exactly 1 bond (k-1)
    for (const auto& kmer : kmers) {
        BOOST_CHECK_EQUAL(kmer.size(), 1);
    }

    std::vector<KmerPath> kmers_k3;
    process_kmers(*mol, 3,
                  [&](const KmerPath& kmer) { kmers_k3.push_back(kmer); });
    // Each k=3 path should have exactly 2 bonds
    for (const auto& kmer : kmers_k3) {
        BOOST_CHECK_EQUAL(kmer.size(), 2);
    }
}

BOOST_AUTO_TEST_CASE(invalid_k_value)
{
    auto mol = helm::helm_to_rdkit("PEPTIDE1{A.C.G}$$$$V2.0");

    BOOST_CHECK_THROW(process_kmers(*mol, 0, [](const KmerPath&) {}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(process_kmers(*mol, 1, [](const KmerPath&) {}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(basic_fingerprints)

BOOST_AUTO_TEST_CASE(simple_generation)
{
    auto fp = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0");

    BOOST_CHECK(fp != nullptr);
    BOOST_CHECK_EQUAL(fp->getNumBits(), 8192); // Default size
    BOOST_CHECK_GT(fp->getNumOnBits(), 0);     // Some bits should be set
}

BOOST_AUTO_TEST_CASE(different_sizes)
{
    BiologicsFingerprintConfig config;

    config.fp_size = 512;
    auto fp512 = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0", config);
    BOOST_CHECK_EQUAL(fp512->getNumBits(), 512);

    config.fp_size = 4096;
    auto fp4096 = generate_fp("PEPTIDE1{A.C.G}$$$$V2.0", config);
    BOOST_CHECK_EQUAL(fp4096->getNumBits(), 4096);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(subset_relationships)

// Test data: (query, target, should_be_subset, description)
const std::vector<std::tuple<std::string, std::string, bool, std::string>>
    subset_test_data = {
        {"PEPTIDE1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Single residue subset"},
        {"PEPTIDE1{T}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Missing residue"},
        {"RNA1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different polymer type"},
        {"CHEM1{A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different polymer type"},
        {"PEPTIDE1{A.C}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Two residue subset (start)"},
        {"PEPTIDE1{C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Two residue subset (middle)"},
        {"PEPTIDE1{A.C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", true,
         "Exact match"},
        {"RNA1{A.C.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Different biologics type"},
        {"PEPTIDE1{A.C.G.A}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Superset (longer)"},
        {"PEPTIDE1{A.G}$$$$V2.0", "PEPTIDE1{A.C.G}$$$$V2.0", false,
         "Non-contiguous (missing C)"},
};

BOOST_DATA_TEST_CASE(subset_detection, bdata::make(subset_test_data), query,
                     target, expected, description)
{
    const bool result = is_substructure_fingerprint_match(query, target);

    BOOST_CHECK_MESSAGE(result == expected, description << ": expected "
                                                        << expected << ", got "
                                                        << result);
}

BOOST_AUTO_TEST_CASE(transitive_subsets)
{
    // Test subset chain: A ⊆ A.C ⊆ A.C.G ⊆ A.C.G.T
    const auto helm_a = "PEPTIDE1{A}$$$$V2.0";
    const auto helm_ac = "PEPTIDE1{A.C}$$$$V2.0";
    const auto helm_acg = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto helm_acgt = "PEPTIDE1{A.C.G.T}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(helm_a, helm_ac));
    BOOST_CHECK(is_substructure_fingerprint_match(helm_ac, helm_acg));
    BOOST_CHECK(is_substructure_fingerprint_match(helm_acg, helm_acgt));
    BOOST_CHECK(
        is_substructure_fingerprint_match(helm_a, helm_acgt)); // Transitive
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(multi_biologics)

BOOST_AUTO_TEST_CASE(disulfide_bond)
{
    // Two chains connected by disulfide
    const auto disulfide =
        "PEPTIDE1{A.C.G}|PEPTIDE2{T.C.E}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$V2.0";
    const auto chain1 = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto chain2 = "PEPTIDE1{T.C.E}$$$$V2.0";

    // Both chains should be represented
    BOOST_CHECK(is_substructure_fingerprint_match(chain1, disulfide));
    BOOST_CHECK(is_substructure_fingerprint_match(chain2, disulfide));
}

BOOST_AUTO_TEST_CASE(cyclic_peptide)
{
    const auto cyclic = "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto linear = "PEPTIDE1{A.C.G.T}$$$$V2.0";
    const auto rotated_cyclic =
        "PEPTIDE1{C.G.T.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto rotated_linear = "PEPTIDE1{C.G.T.A}$$$$V2.0";

    auto fp_cyclic = generate_fp(cyclic);
    auto fp_linear = generate_fp(linear);
    auto fp_rotated_cyclic = generate_fp(rotated_cyclic);
    auto fp_rotated_linear = generate_fp(rotated_linear);

    // Cyclic should have at least as many features as linear
    BOOST_CHECK_GE(fp_cyclic->getNumOnBits(), fp_linear->getNumOnBits());
    BOOST_CHECK_GE(fp_rotated_cyclic->getNumOnBits(),
                   fp_rotated_linear->getNumOnBits());

    // all variations  should be subset of cyclic
    BOOST_CHECK(is_substructure_fingerprint_match(linear, cyclic));
    BOOST_CHECK(is_substructure_fingerprint_match(rotated_cyclic, cyclic));
    BOOST_CHECK(is_substructure_fingerprint_match(rotated_linear, cyclic));

    BOOST_CHECK(is_substructure_fingerprint_match(linear, rotated_cyclic));
    BOOST_CHECK(is_substructure_fingerprint_match(cyclic, rotated_cyclic));
    BOOST_CHECK(
        is_substructure_fingerprint_match(rotated_linear, rotated_cyclic));
}

BOOST_AUTO_TEST_CASE(peptide_rna_conjugate)
{
    const auto conjugate =
        "PEPTIDE1{A.C.G}|RNA1{A.C.G}$PEPTIDE1,RNA1,3:R2-1:R1$$$V2.0";
    const auto peptide = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto rna = "RNA1{A.C.G}$$$$V2.0";

    // Both biologics types should be represented
    BOOST_CHECK(is_substructure_fingerprint_match(peptide, conjugate));
    BOOST_CHECK(is_substructure_fingerprint_match(rna, conjugate));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(repeat_sequences)

BOOST_AUTO_TEST_CASE(homopolymer_limitation)
{
    // CRITICAL LIMITATION: Sequences longer than max_k with same residue
    // have identical fingerprints.

    BiologicsFingerprintConfig config;
    config.max_k = 4;

    auto fp_4 = generate_fp("PEPTIDE1{A.A.A.A}$$$$V2.0", config);
    auto fp_5 = generate_fp("PEPTIDE1{A.A.A.A.A}$$$$V2.0", config);
    auto fp_10 = generate_fp("PEPTIDE1{A.A.A.A.A.A.A.A.A.A}$$$$V2.0", config);

    // All should have identical k-mer sets when k ≤ max_k
    // This is a known limitation, not a bug
    BOOST_CHECK_EQUAL(fp_4->getNumOnBits(), fp_5->getNumOnBits());
    BOOST_CHECK_EQUAL(fp_4->getNumOnBits(), fp_10->getNumOnBits());
}

BOOST_AUTO_TEST_CASE(context_helps_distinguish_repeats)
{
    // Flanking sequences create unique k-mers that help distinguish
    // different repeat counts (up to max_k)

    auto helm_2a = "PEPTIDE1{T.A.A.E}$$$$V2.0";
    auto helm_3a = "PEPTIDE1{T.A.A.A.E}$$$$V2.0";

    // helm_2a should not be subset of helm_3a (fewer As)
    BOOST_CHECK(!is_substructure_fingerprint_match(helm_2a, helm_3a));

    // helm_3a should NOT be subset of helm_2a (has unique T-A-A-A 4-mer)
    BOOST_CHECK(!is_substructure_fingerprint_match(helm_3a, helm_2a));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(branching_structures)

BOOST_AUTO_TEST_CASE(simple_single_branch)
{
    // Single branch point: A.C(G)T
    // Note: HELM allows only one branch monomer
    const auto branched = "PEPTIDE1{A.C(G)T}$$$$V2.0";
    const auto act = "PEPTIDE1{A.C.T}$$$$V2.0";
    const auto acg_linear = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto acg_branched = "PEPTIDE1{A.C(G)}$$$$V2.0";

    auto fp_branched = generate_fp(branched);
    auto fp_act = generate_fp(act);

    // Linear paths should be subsets
    BOOST_CHECK(is_substructure_fingerprint_match(act, branched));
    BOOST_CHECK(!is_substructure_fingerprint_match(acg_linear, branched));
    BOOST_CHECK(is_substructure_fingerprint_match(acg_branched, branched));

    // Branched structure has both paths
    BOOST_CHECK_GT(fp_branched->getNumOnBits(), fp_act->getNumOnBits());
}

BOOST_AUTO_TEST_CASE(multi_polymer_y_branch)
{
    // Y-shaped branch using multi-polymer: A.C branches to G.K and E.F
    // PEPTIDE1{A(C.G.K)P} becomes two polymers with custom connection
    const auto y_branch =
        "PEPTIDE1{A.C}|PEPTIDE2{G.K}|PEPTIDE3{E.F}$"
        "PEPTIDE1,PEPTIDE2,2:R2-1:R1|PEPTIDE1,PEPTIDE3,2:R3-1:R1$$$V2.0";

    // All three arms should be represented
    const auto arm1 = "PEPTIDE1{A.C}$$$$V2.0";
    const auto arm2 = "PEPTIDE1{G.K}$$$$V2.0";
    const auto arm3 = "PEPTIDE1{E.F}$$$$V2.0";

    const auto joined1 = "PEPTIDE1{A.C.G}$$$$V2.0";
    const auto joined2 = "PEPTIDE1{C.G.K}$$$$V2.0";
    const auto joined3 = "PEPTIDE1{A.C.G.K}$$$$V2.0";
    const auto joined4 = "PEPTIDE1{A.C(E)G.K}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(arm1, y_branch));
    BOOST_CHECK(is_substructure_fingerprint_match(arm2, y_branch));
    BOOST_CHECK(is_substructure_fingerprint_match(arm3, y_branch));

    BOOST_CHECK(is_substructure_fingerprint_match(joined1, y_branch));
    BOOST_CHECK(is_substructure_fingerprint_match(joined2, y_branch));
    BOOST_CHECK(is_substructure_fingerprint_match(joined3, y_branch));
    BOOST_CHECK(is_substructure_fingerprint_match(joined4, y_branch));
}

BOOST_AUTO_TEST_CASE(multi_polymer_star_topology)
{
    // Star: central hub K with three arms
    const auto star = "PEPTIDE1{K}|PEPTIDE2{A.C}|PEPTIDE3{G.T}|PEPTIDE4{E.F}$"
                      "PEPTIDE1,PEPTIDE2,1:R1-1:R1|"
                      "PEPTIDE1,PEPTIDE3,1:R2-1:R1|"
                      "PEPTIDE1,PEPTIDE4,1:R3-1:R1$$$V2.0";

    // Simple arms
    const auto ac = "PEPTIDE1{A.C}$$$$V2.0";
    const auto gt = "PEPTIDE1{G.T}$$$$V2.0";
    const auto ef = "PEPTIDE1{E.F}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(ac, star));
    BOOST_CHECK(is_substructure_fingerprint_match(gt, star));
    BOOST_CHECK(is_substructure_fingerprint_match(ef, star));

    // Complex: paths through the hub using multi-polymer
    const auto path_ac_gt =
        "PEPTIDE1{K}|PEPTIDE2{A.C}|PEPTIDE3{G.T}$"
        "PEPTIDE1,PEPTIDE2,1:R1-1:R1|PEPTIDE1,PEPTIDE3,1:R2-1:R1$$$V2.0";
    const auto path_all =
        "PEPTIDE1{K}|PEPTIDE2{A.C}|PEPTIDE3{G.T}|PEPTIDE4{E.F}$"
        "PEPTIDE1,PEPTIDE2,1:R1-1:R1|PEPTIDE1,PEPTIDE3,1:R2-1:R1|"
        "PEPTIDE1,PEPTIDE4,1:R3-1:R1$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(path_ac_gt, star));
    BOOST_CHECK(is_substructure_fingerprint_match(path_all, star));
}

BOOST_AUTO_TEST_CASE(ladder_structure)
{
    // Ladder: two parallel chains with cross-links at positions 2 and 3
    const auto ladder =
        "PEPTIDE1{A.C.G.T}|PEPTIDE2{E.F.H.I}$"
        "PEPTIDE1,PEPTIDE2,2:R3-2:R3|PEPTIDE1,PEPTIDE2,3:R3-3:R3$$$V2.0";

    // Simple rails
    const auto rail1 = "PEPTIDE1{A.C.G.T}$$$$V2.0";
    const auto rail2 = "PEPTIDE1{E.F.H.I}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(rail1, ladder));
    BOOST_CHECK(is_substructure_fingerprint_match(rail2, ladder));

    // Complex: partial ladder with one cross-link
    const auto partial_ladder = "PEPTIDE1{A.C.G}|PEPTIDE2{E.F.H}$"
                                "PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial_ladder, ladder));

    // Complex: sub-ladder with both cross-links
    const auto sub_ladder =
        "PEPTIDE1{C.G.T}|PEPTIDE2{F.H.I}$"
        "PEPTIDE1,PEPTIDE2,1:R3-1:R3|PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(sub_ladder, ladder));

    // Negative: wrong cross-link positions should not match
    const auto wrong_links =
        "PEPTIDE1{A.C.G.T}|PEPTIDE2{E.F.H.I}$"
        "PEPTIDE1,PEPTIDE2,1:R3-1:R3|PEPTIDE1,PEPTIDE2,4:R3-4:R3$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_links, ladder));
}

BOOST_AUTO_TEST_CASE(tree_hierarchy)
{
    // Tree: A branches to C and E, where E further branches to G and K
    const auto tree = "PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{E}|"
                      "PEPTIDE4{G}|PEPTIDE5{K}$"
                      "PEPTIDE1,PEPTIDE2,1:R1-1:R1|"
                      "PEPTIDE1,PEPTIDE3,1:R2-1:R1|"
                      "PEPTIDE3,PEPTIDE4,1:R1-1:R1|"
                      "PEPTIDE3,PEPTIDE5,1:R2-1:R1$$$V2.0";

    // Simple: individual branches
    const auto path_ac =
        "PEPTIDE1{A}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:R1-1:R1$$$V2.0";
    const auto path_ae =
        "PEPTIDE1{A}|PEPTIDE2{E}$PEPTIDE1,PEPTIDE2,1:R2-1:R1$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(path_ac, tree));
    BOOST_CHECK(is_substructure_fingerprint_match(path_ae, tree));

    // Complex: paths through multiple levels (A->E->G)
    const auto path_aeg =
        "PEPTIDE1{A}|PEPTIDE2{E}|PEPTIDE3{G}$"
        "PEPTIDE1,PEPTIDE2,1:R2-1:R1|PEPTIDE2,PEPTIDE3,1:R1-1:R1$$$V2.0";
    const auto path_aek =
        "PEPTIDE1{A}|PEPTIDE2{E}|PEPTIDE3{K}$"
        "PEPTIDE1,PEPTIDE2,1:R2-1:R1|PEPTIDE2,PEPTIDE3,1:R2-1:R1$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(path_aeg, tree));
    BOOST_CHECK(is_substructure_fingerprint_match(path_aek, tree));

    // Complex: subtree with all E branches
    const auto subtree_e =
        "PEPTIDE1{E}|PEPTIDE2{G}|PEPTIDE3{K}$"
        "PEPTIDE1,PEPTIDE2,1:R1-1:R1|PEPTIDE1,PEPTIDE3,1:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(subtree_e, tree));

    // Negative: wrong connections
    const auto wrong_tree =
        "PEPTIDE1{A}|PEPTIDE2{G}$PEPTIDE1,PEPTIDE2,1:R1-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_tree, tree));
}

BOOST_AUTO_TEST_CASE(multiple_disulfide_bonds)
{
    // Insulin-like: two chains with multiple disulfide bridges
    const auto insulin_like =
        "PEPTIDE1{G.I.V.E.C.C.A}|PEPTIDE2{F.V.N.C.H.L.C.G}$"
        "PEPTIDE1,PEPTIDE2,5:R3-4:R3|PEPTIDE1,PEPTIDE2,6:R3-7:R3$$$V2.0";

    // Simple: individual chains
    const auto chain1 = "PEPTIDE1{G.I.V.E.C.C.A}$$$$V2.0";
    const auto chain2 = "PEPTIDE1{F.V.N.C.H.L.C.G}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(chain1, insulin_like));
    BOOST_CHECK(is_substructure_fingerprint_match(chain2, insulin_like));

    // Complex: subsequences from both chains
    const auto sub1 = "PEPTIDE1{I.V.E.C}$$$$V2.0";
    const auto sub2 = "PEPTIDE1{N.C.H.L}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(sub1, insulin_like));
    BOOST_CHECK(is_substructure_fingerprint_match(sub2, insulin_like));

    // Complex: structure with both chains but only one disulfide
    const auto single_ss = "PEPTIDE1{G.I.V.E.C.C.A}|PEPTIDE2{F.V.N.C.H.L.C.G}$"
                           "PEPTIDE1,PEPTIDE2,5:R3-4:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(single_ss, insulin_like));

    // Complex: partial chains with disulfide
    const auto partial =
        "PEPTIDE1{E.C.C}|PEPTIDE2{C.H.L.C}$"
        "PEPTIDE1,PEPTIDE2,2:R3-1:R3|PEPTIDE1,PEPTIDE2,3:R3-4:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial, insulin_like));

    // Negative: wrong disulfide positions
    const auto wrong_ss = "PEPTIDE1{G.I.V.E.C.C.A}|PEPTIDE2{F.V.N.C.H.L.C.G}$"
                          "PEPTIDE1,PEPTIDE2,5:R3-7:R3$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_ss, insulin_like));
}

BOOST_AUTO_TEST_CASE(multi_type_conjugate_chain)
{
    // Peptide-RNA-Peptide chain
    const auto conjugate =
        "PEPTIDE1{A.C}|RNA1{G.U}|PEPTIDE2{T.E}$"
        "PEPTIDE1,RNA1,2:R2-1:R1|RNA1,PEPTIDE2,2:R2-1:R1$$$V2.0";

    // Simple: individual segments
    const auto pep1 = "PEPTIDE1{A.C}$$$$V2.0";
    const auto rna = "RNA1{G.U}$$$$V2.0";
    const auto pep2 = "PEPTIDE1{T.E}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(pep1, conjugate));
    BOOST_CHECK(is_substructure_fingerprint_match(rna, conjugate));
    BOOST_CHECK(is_substructure_fingerprint_match(pep2, conjugate));

    // Complex: two-segment paths
    const auto pep_rna =
        "PEPTIDE1{A.C}|RNA1{G.U}$PEPTIDE1,RNA1,2:R2-1:R1$$$V2.0";
    const auto rna_pep =
        "RNA1{G.U}|PEPTIDE1{T.E}$RNA1,PEPTIDE1,2:R2-1:R1$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(pep_rna, conjugate));
    BOOST_CHECK(is_substructure_fingerprint_match(rna_pep, conjugate));

    // Complex: partial three-segment
    const auto partial =
        "PEPTIDE1{C}|RNA1{G.U}|PEPTIDE2{T}$"
        "PEPTIDE1,RNA1,1:R2-1:R1|RNA1,PEPTIDE2,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial, conjugate));

    // Negative: wrong order
    const auto wrong_order =
        "RNA1{G.U}|PEPTIDE1{A.C}|PEPTIDE2{T.E}$"
        "RNA1,PEPTIDE1,1:R1-1:R1|PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_order, conjugate));
}

BOOST_AUTO_TEST_CASE(circular_with_branch)
{
    // Cyclic peptide with a branch off the ring at position 2 (C)
    const auto cyclic_branched =
        "PEPTIDE1{A.C.G.T}|PEPTIDE2{K.E}$"
        "PEPTIDE1,PEPTIDE1,4:R2-1:R1|PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0";

    // Simple: backbone and branch separately
    const auto cyclic = "PEPTIDE1{A.C.G.T}$$$$V2.0";
    const auto branch = "PEPTIDE1{K.E}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(cyclic, cyclic_branched));
    BOOST_CHECK(is_substructure_fingerprint_match(branch, cyclic_branched));

    // Complex: partial ring with branch point
    const auto partial_ring_branch =
        "PEPTIDE1{A.C.G}|PEPTIDE2{K.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial_ring_branch,
                                                  cyclic_branched));

    // Complex: cyclic backbone with branch attachment point
    const auto ring_with_attach =
        "PEPTIDE1{A.C.G.T}|PEPTIDE2{K}$"
        "PEPTIDE1,PEPTIDE1,4:R2-1:R1|PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0";
    BOOST_CHECK(
        is_substructure_fingerprint_match(ring_with_attach, cyclic_branched));

    // Complex: with branching annotation
    const auto with_notation = "PEPTIDE1{A.C(K)G.T}$$$$V2.0";
    BOOST_CHECK(
        is_substructure_fingerprint_match(with_notation, cyclic_branched));

    // Negative: branch at wrong position
    const auto wrong_pos =
        "PEPTIDE1{A.C.G.T}|PEPTIDE2{K.E}$"
        "PEPTIDE1,PEPTIDE1,4:R2-1:R1|PEPTIDE1,PEPTIDE2,1:R3-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_pos, cyclic_branched));
}

BOOST_AUTO_TEST_CASE(multiple_connection_points)
{
    // Two chains with connections at positions 1, 3, and 5
    const auto multi_connect = "PEPTIDE1{A.C.G.T.E}|PEPTIDE2{K.F.H.I.L}$"
                               "PEPTIDE1,PEPTIDE2,1:R3-1:R3|"
                               "PEPTIDE1,PEPTIDE2,3:R3-3:R3|"
                               "PEPTIDE1,PEPTIDE2,5:R3-5:R3$$$V2.0";

    // Simple: individual chains
    const auto chain1 = "PEPTIDE1{A.C.G.T.E}$$$$V2.0";
    const auto chain2 = "PEPTIDE1{K.F.H.I.L}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(chain1, multi_connect));
    BOOST_CHECK(is_substructure_fingerprint_match(chain2, multi_connect));

    // Complex: with subset of connections
    const auto two_connect =
        "PEPTIDE1{A.C.G.T.E}|PEPTIDE2{K.F.H.I.L}$"
        "PEPTIDE1,PEPTIDE2,1:R3-1:R3|PEPTIDE1,PEPTIDE2,3:R3-3:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(two_connect, multi_connect));

    // Complex: single connection
    const auto single_connect = "PEPTIDE1{A.C.G.T.E}|PEPTIDE2{K.F.H.I.L}$"
                                "PEPTIDE1,PEPTIDE2,1:R3-1:R3$$$V2.0";
    BOOST_CHECK(
        is_substructure_fingerprint_match(single_connect, multi_connect));

    // Complex: partial chains with connections
    const auto partial =
        "PEPTIDE1{A.C.G}|PEPTIDE2{K.F.H}$"
        "PEPTIDE1,PEPTIDE2,1:R3-1:R3|PEPTIDE1,PEPTIDE2,3:R3-3:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial, multi_connect));

    // Verify different connection counts produce different fingerprints
    auto fp_multi = generate_fp(multi_connect);
    auto fp_single_connect = generate_fp(single_connect);
    BOOST_CHECK_NE(fp_multi->getNumOnBits(), fp_single_connect->getNumOnBits());

    // Negative: wrong connection positions
    const auto wrong_connect =
        "PEPTIDE1{A.C.G.T.E}|PEPTIDE2{K.F.H.I.L}$"
        "PEPTIDE1,PEPTIDE2,2:R3-2:R3|PEPTIDE1,PEPTIDE2,4:R3-4:R3$$$V2.0";
    BOOST_CHECK(
        !is_substructure_fingerprint_match(wrong_connect, multi_connect));
}

BOOST_AUTO_TEST_CASE(symmetric_dumbbell)
{
    // Symmetric: two identical A.C.G arms connected at their ends
    const auto dumbbell = "PEPTIDE1{A.C.G}|PEPTIDE2{A.C.G}$"
                          "PEPTIDE1,PEPTIDE2,3:R2-1:R1$$$V2.0";
    const auto arm = "PEPTIDE1{A.C.G}$$$$V2.0";

    auto fp_dumbbell = generate_fp(dumbbell);
    auto fp_arm = generate_fp(arm);

    // Simple: single arm
    BOOST_CHECK(is_substructure_fingerprint_match(arm, dumbbell));

    // Complex: partial dumbbell (one complete arm + connection point)
    const auto partial_dumbbell =
        "PEPTIDE1{A.C.G}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,3:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial_dumbbell, dumbbell));

    // Complex: both arms with partial sequences
    const auto both_partial =
        "PEPTIDE1{C.G}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(both_partial, dumbbell));

    // Negative: single branched arm with same sequence
    const auto branched_arm = "PEPTIDE1{A.C(G)}$$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(branched_arm, dumbbell));

    // Verify connection adds features
    BOOST_CHECK_GE(fp_dumbbell->getNumOnBits(), fp_arm->getNumOnBits());

    // Negative: asymmetric structure
    const auto asymmetric =
        "PEPTIDE1{A.C.G}|PEPTIDE2{T.E.F}$PEPTIDE1,PEPTIDE2,3:R2-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(asymmetric, dumbbell));

    // Negative: wrong connection point
    const auto wrong_connection =
        "PEPTIDE1{A.C.G}|PEPTIDE2{A.C.G}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_connection, dumbbell));
}

BOOST_AUTO_TEST_CASE(cross_polymer_bridge)
{
    // Bridge: Peptide-RNA-Peptide with specific connections
    const auto bridge =
        "PEPTIDE1{A.C}|RNA1{G.U}|PEPTIDE2{T.E}$"
        "PEPTIDE1,RNA1,2:R3-1:R3|RNA1,PEPTIDE2,2:R3-1:R3$$$V2.0";

    // Simple: individual components
    const auto pep1 = "PEPTIDE1{A.C}$$$$V2.0";
    const auto rna = "RNA1{G.U}$$$$V2.0";
    const auto pep2 = "PEPTIDE1{T.E}$$$$V2.0";

    auto fp_bridge = generate_fp(bridge);
    auto fp_pep1 = generate_fp(pep1);
    auto fp_rna = generate_fp(rna);

    BOOST_CHECK(is_substructure_fingerprint_match(pep1, bridge));
    BOOST_CHECK(is_substructure_fingerprint_match(rna, bridge));
    BOOST_CHECK(is_substructure_fingerprint_match(pep2, bridge));

    // Complex: two-component bridges
    const auto pep_rna_bridge =
        "PEPTIDE1{A.C}|RNA1{G.U}$PEPTIDE1,RNA1,2:R3-1:R3$$$V2.0";
    const auto rna_pep_bridge =
        "RNA1{G.U}|PEPTIDE1{T.E}$RNA1,PEPTIDE1,2:R3-1:R3$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(pep_rna_bridge, bridge));
    BOOST_CHECK(is_substructure_fingerprint_match(rna_pep_bridge, bridge));

    // Complex: partial sequences with cross-type connection
    const auto partial_cross =
        "PEPTIDE1{C}|RNA1{G}$PEPTIDE1,RNA1,1:R3-1:R3$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial_cross, bridge));

    // Verify cross-type features present
    BOOST_CHECK_GT(fp_bridge->getNumOnBits(), fp_pep1->getNumOnBits());
    BOOST_CHECK_GT(fp_bridge->getNumOnBits(), fp_rna->getNumOnBits());

    // Negative: different connection pattern
    const auto wrong_bridge =
        "PEPTIDE1{A.C}|RNA1{G.U}|PEPTIDE2{T.E}$"
        "PEPTIDE1,RNA1,1:R3-1:R3|RNA1,PEPTIDE2,1:R3-1:R3$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_bridge, bridge));
}

BOOST_AUTO_TEST_CASE(cage_structure)
{
    // Cage: three chains forming a triangular loop
    // A.C connects to G.T, G.T connects to E.F, E.F connects back to A.C
    const auto cage = "PEPTIDE1{A.C}|PEPTIDE2{G.T}|PEPTIDE3{E.F}$"
                      "PEPTIDE1,PEPTIDE2,2:R2-1:R1|"
                      "PEPTIDE2,PEPTIDE3,2:R2-1:R1|"
                      "PEPTIDE3,PEPTIDE1,2:R2-1:R1$$$V2.0";

    // Simple: individual edges
    const auto edge1 = "PEPTIDE1{A.C}$$$$V2.0";
    const auto edge2 = "PEPTIDE1{G.T}$$$$V2.0";
    const auto edge3 = "PEPTIDE1{E.F}$$$$V2.0";

    BOOST_CHECK(is_substructure_fingerprint_match(edge1, cage));
    BOOST_CHECK(is_substructure_fingerprint_match(edge2, cage));
    BOOST_CHECK(is_substructure_fingerprint_match(edge3, cage));

    // Complex: two connected edges (partial cage)
    const auto two_edges =
        "PEPTIDE1{A.C}|PEPTIDE2{G.T}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(two_edges, cage));

    // Complex: all three edges but only two connections (not closed)
    const auto open_cage =
        "PEPTIDE1{A.C}|PEPTIDE2{G.T}|PEPTIDE3{E.F}$"
        "PEPTIDE1,PEPTIDE2,2:R2-1:R1|PEPTIDE2,PEPTIDE3,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(open_cage, cage));

    // Complex: partial sequences with connections
    const auto partial_cage =
        "PEPTIDE1{C}|PEPTIDE2{G.T}|PEPTIDE3{E}$"
        "PEPTIDE1,PEPTIDE2,1:R2-1:R1|PEPTIDE2,PEPTIDE3,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(is_substructure_fingerprint_match(partial_cage, cage));

    // Negative: wrong connection topology
    const auto wrong_topology =
        "PEPTIDE1{A.C}|PEPTIDE2{G.T}|PEPTIDE3{E.F}$"
        "PEPTIDE1,PEPTIDE2,2:R2-1:R1|PEPTIDE1,PEPTIDE3,2:R2-1:R1$$$V2.0";
    BOOST_CHECK(!is_substructure_fingerprint_match(wrong_topology, cage));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(cycle_features)

BOOST_AUTO_TEST_CASE(linear_vs_cyclic_fingerprints)
{
    const auto linear = "PEPTIDE1{A.C.G.T}$$$$V2.0";
    const auto cyclic = "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";

    auto fp_linear = generate_fp(linear);
    auto fp_cyclic = generate_fp(cyclic);

    // Cyclic and linear fingerprints should be different
    BOOST_CHECK(*fp_linear != *fp_cyclic);

    // Linear should be subset of cyclic (all linear k-mers are present)
    BOOST_CHECK(is_substructure_fingerprint_match(linear, cyclic));

    // Cyclic should NOT be subset of linear (has unique cycle features)
    BOOST_CHECK(!is_substructure_fingerprint_match(cyclic, linear));
}

BOOST_AUTO_TEST_CASE(different_cycle_sizes)
{
    // Verify cycles of different sizes produce distinct fingerprints
    const auto cycle4 = "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle6 =
        "PEPTIDE1{A.C.G.T.E.F}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0";
    const auto cycle8 = "PEPTIDE1{A.C.G.T.E.F.H.I}$"
                        "PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0";

    auto fp_cycle4 = generate_fp(cycle4);
    auto fp_cycle6 = generate_fp(cycle6);
    auto fp_cycle8 = generate_fp(cycle8);

    // Different cycle sizes should have different fingerprints
    BOOST_CHECK(*fp_cycle4 != *fp_cycle6);
    BOOST_CHECK(*fp_cycle4 != *fp_cycle8);
    BOOST_CHECK(*fp_cycle6 != *fp_cycle8);
}

BOOST_AUTO_TEST_CASE(cycle_composition_matters)
{
    // Cycles with different residue composition should produce different
    // fingerprints
    const auto cycle_acgt =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_aagg =
        "PEPTIDE1{A.A.G.G}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_eeee =
        "PEPTIDE1{E.E.E.E}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";

    auto fp_acgt = generate_fp(cycle_acgt);
    auto fp_aagg = generate_fp(cycle_aagg);
    auto fp_eeee = generate_fp(cycle_eeee);

    // Different compositions should produce different fingerprints
    BOOST_CHECK(*fp_acgt != *fp_aagg);
    BOOST_CHECK(*fp_acgt != *fp_eeee);
    BOOST_CHECK(*fp_aagg != *fp_eeee);
}

BOOST_AUTO_TEST_CASE(rotation_invariance)
{
    // Rotated representations of the same cycle should have identical
    // fingerprints
    const auto cycle_acgt =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_cgta =
        "PEPTIDE1{C.G.T.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_gtac =
        "PEPTIDE1{G.T.A.C}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_tacg =
        "PEPTIDE1{T.A.C.G}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";

    auto fp_acgt = generate_fp(cycle_acgt);
    auto fp_cgta = generate_fp(cycle_cgta);
    auto fp_gtac = generate_fp(cycle_gtac);
    auto fp_tacg = generate_fp(cycle_tacg);

    // All rotations should have identical fingerprints
    BOOST_CHECK(*fp_acgt == *fp_cgta);
    BOOST_CHECK(*fp_acgt == *fp_gtac);
    BOOST_CHECK(*fp_acgt == *fp_tacg);
}

BOOST_AUTO_TEST_CASE(multiple_cycles)
{
    // Structure with two separate cycles should capture both
    const auto two_cycles =
        "PEPTIDE1{A.C.G}|PEPTIDE2{T.E.F}$"
        "PEPTIDE1,PEPTIDE1,3:R2-1:R1|PEPTIDE2,PEPTIDE2,3:R2-1:R1$$$V2.0";

    const auto single_cycle1 =
        "PEPTIDE1{A.C.G}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0";
    const auto single_cycle2 =
        "PEPTIDE1{T.E.F}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0";

    // Both individual cycles should be substructures of the dual-cycle
    // structure
    BOOST_CHECK(is_substructure_fingerprint_match(single_cycle1, two_cycles));
    BOOST_CHECK(is_substructure_fingerprint_match(single_cycle2, two_cycles));
}

BOOST_AUTO_TEST_CASE(different_linkage_types)
{
    // Same residues, same cycle size, but different attachment points
    const auto cycle_r2_r1 =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_r3_r1 =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R3-1:R1$$$V2.0";
    const auto cycle_r3_r3 =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R3-1:R3$$$V2.0";

    auto fp_r2_r1 = generate_fp(cycle_r2_r1);
    auto fp_r3_r1 = generate_fp(cycle_r3_r1);
    auto fp_r3_r3 = generate_fp(cycle_r3_r3);

    // Different linkage types should NOT match as substructures
    BOOST_CHECK(!is_substructure_fingerprint_match(cycle_r2_r1, cycle_r3_r1));
    BOOST_CHECK(!is_substructure_fingerprint_match(cycle_r3_r1, cycle_r2_r1));
    BOOST_CHECK(!is_substructure_fingerprint_match(cycle_r2_r1, cycle_r3_r3));
}

BOOST_AUTO_TEST_CASE(linkage_with_same_linear_sequence)
{
    // Verify that cycles with different linkages don't match even if their
    // linear sequence is identical
    const auto linear = "PEPTIDE1{A.C.G.T}$$$$V2.0";
    const auto cycle_r2_r1 =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    const auto cycle_r3_r1 =
        "PEPTIDE1{A.C.G.T}$PEPTIDE1,PEPTIDE1,4:R3-1:R1$$$V2.0";

    // Linear should be subset of both cycles (same k-mers)
    BOOST_CHECK(is_substructure_fingerprint_match(linear, cycle_r2_r1));
    BOOST_CHECK(is_substructure_fingerprint_match(linear, cycle_r3_r1));

    // But the two cycles should NOT match each other (different cyclization)
    BOOST_CHECK(!is_substructure_fingerprint_match(cycle_r2_r1, cycle_r3_r1));
    BOOST_CHECK(!is_substructure_fingerprint_match(cycle_r3_r1, cycle_r2_r1));
}

BOOST_AUTO_TEST_CASE(canonical_rotation_with_repeats)
{
    // Test edge case with repeated elements
    auto cyclic_aabc = "PEPTIDE1{A.A.B.C}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    auto cyclic_abca = "PEPTIDE1{A.B.C.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    auto cyclic_bcaa = "PEPTIDE1{B.C.A.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";
    auto cyclic_caab = "PEPTIDE1{C.A.A.B}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0";

    auto fp1 = generate_fp(cyclic_aabc);
    auto fp2 = generate_fp(cyclic_abca);
    auto fp3 = generate_fp(cyclic_bcaa);
    auto fp4 = generate_fp(cyclic_caab);

    // All should be equivalent
    BOOST_CHECK(*fp1 == *fp2);
    BOOST_CHECK(*fp2 == *fp3);
    BOOST_CHECK(*fp3 == *fp4);
}

BOOST_AUTO_TEST_SUITE_END()
