
#define BOOST_TEST_MODULE monomeric

#include <algorithm>
#include <array>
#include <unordered_set>

#include <rdkit/GraphMol/RWMol.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <fmt/format.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

using namespace boost::unit_test;

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_validate_monomers)
{
    for (const auto& helm :
         {"PEPTIDE1{A.C.W}$$$$V2.0", "RNA1{R(A)P.[dR](C)P}$$$$V2.0",
          "PEPTIDE1{X}$$$$V2.0", "PEPTIDE1{A.X.W}$$$$V2.0",
          "RNA1{R(N)P.[dR](N)P}$$$$V2.0", "CHEM1{[CCO]}$$$$V2.0",
          "PEPTIDE1{[C* |$;_R1$|]}$$$$V2.0"}) {
        auto mol = rdkit_extensions::to_rdkit(helm);
        BOOST_CHECK_NO_THROW(validate_monomers(*mol));
    }

    for (const auto& [helm, expected] :
         std::vector<std::pair<std::string, std::string>>{
             {"PEPTIDE1{A.[missingMonomer]}$$$$V2.0",
              "Peptide monomer missingMonomer not found in monomer database"},
             {"RNA1{R([missingMonomer])P}$$$$V2.0",
              "Nucleic acid monomer missingMonomer not found in monomer "
              "database"},
             {"CHEM1{[missingMonomer]}$$$$V2.0",
              "CHEM monomer missingMonomer not found in monomer database"},
             {"CHEM1{W}$$$$V2.0",
              "CHEM monomer W not found in monomer database"},
             {"CHEM1{X}$$$$V2.0",
              "CHEM monomer X not found in monomer database"},
             {"CHEM1{N}$$$$V2.0",
              "CHEM monomer N not found in monomer database"},
             {"RNA1{R(X)P}$$$$V2.0",
              "Nucleic acid monomer X not found in monomer database"},
             {"PEPTIDE1{X.[missingMonomer]}$$$$V2.0",
              "Peptide monomer missingMonomer not found in monomer database"},
             {"RNA1{R(N)P.R([missingMonomer])P}$$$$V2.0",
              "Nucleic acid monomer missingMonomer not found in monomer "
              "database"}}) {
        auto mol = rdkit_extensions::to_rdkit(helm);
        auto check_error = [&](const std::runtime_error& error) {
            return error.what() == expected;
        };
        BOOST_CHECK_EXCEPTION(validate_monomers(*mol), std::runtime_error,
                              check_error);
        // Missing SMILES_MONOMER properties also mean database monomers.
        for (auto* monomer : mol->atoms()) {
            monomer->clearProp(SMILES_MONOMER);
        }
        BOOST_CHECK_EXCEPTION(validate_monomers(*mol), std::runtime_error,
                              check_error);
    }

    auto mol = rdkit_extensions::to_rdkit("CHEM1{[CCO]}$$$$V2.0");
    auto* monomer = mol->getAtomWithIdx(0);
    monomer->setProp(ATOM_LABEL, std::string("C1CC"));
    BOOST_CHECK_EXCEPTION(validate_monomers(*mol), std::runtime_error,
                          [](const std::runtime_error& error) {
                              return std::string(error.what()) ==
                                     "Could not parse monomer SMILES: C1CC";
                          });
}

BOOST_AUTO_TEST_CASE(test_get_attachment_points_for_smiles)
{
    const std::array<std::string, 5> alanine_smiles = {
        "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB "
        ":1.pdbName. CA :2.pdbName. N  :3.pdbName. H  :4.pdbName. C  "
        ":5.pdbName. O  :6.pdbName. OXT|",
        "C[C@H](N[H:1])C(=O)[OH:2]", "[1*]N[C@@H](C)C(=O)O[2*]",
        "*N[C@@H](C)C(=O)O* |$_R1;;;;;;;_R2$,atomProp:1.pdbName. N  "
        ":2.pdbName. CA :3.pdbName. CB :4.pdbName. C  :5.pdbName. O  "
        ",a:2|",
        "*N[C@@H](C)C(=O)O* |$_R1;;;;;;;_R2$|"};
    const std::vector<std::pair<int, std::string>> expected = {{1, "N"},
                                                               {2, "O"}};

    for (const auto& smiles : alanine_smiles) {
        BOOST_TEST(get_attachment_points_for_smiles(smiles) == expected);
    }
}

BOOST_AUTO_TEST_CASE(test_normalize_smiles_attachment_points)
{
    const std::string alanine = "*N[C@@H](C)C(=O)* |$_R1;;;;;;_R2$|";
    const std::vector<std::tuple<std::string, std::string,
                                 std::vector<std::pair<int, std::string>>>>
        test_cases = {
            {"C[C@H](N[H:1])C(=O)[OH:2]", alanine, {{1, "N"}, {2, "C"}}},
            {"[1*]N[C@@H](C)C(=O)[2*]", alanine, {{1, "N"}, {2, "C"}}},
            {alanine, alanine, {{1, "N"}, {2, "C"}}},
            {"[*:1]N[C@@H](C)C(=O)[*:2]", alanine, {{1, "N"}, {2, "C"}}},
            {"C[C@H](N)C(=O)O |$;;_R1;;;_R2$|",
             "*N[C@@H](C)C(=O)O* |$_R1;;;;;;;_R2$|",
             {{1, "N"}, {2, "O"}}},
            {"O=P(O)([OH:1])[OH:2]",
             "O=P(O)(*)* |$;;;_R1;_R2$|",
             {{1, "P"}, {2, "P"}}}};

    for (const auto& [input_smiles, expected_smiles, expected_aps] :
         test_cases) {
        const auto normalized =
            normalize_smiles_attachment_points(input_smiles);
        BOOST_TEST(normalized.find("_R1") != std::string::npos);
        BOOST_TEST(normalized.find("_R2") != std::string::npos);
        BOOST_TEST(get_attachment_points_for_smiles(normalized) ==
                   expected_aps);
        const auto mol = rdkit_extensions::to_rdkit(
            normalized, rdkit_extensions::Format::EXTENDED_SMILES);
        const auto expected = rdkit_extensions::to_rdkit(
            expected_smiles, rdkit_extensions::Format::EXTENDED_SMILES);
        BOOST_TEST(rdkit_extensions::to_string(
                       *mol, rdkit_extensions::Format::EXTENDED_SMILES) ==
                   rdkit_extensions::to_string(
                       *expected, rdkit_extensions::Format::EXTENDED_SMILES));

        std::vector<unsigned int> attachment_point_nums;
        for (const auto* atom : mol->atoms()) {
            const auto r_group_num = rdkit_extensions::get_r_group_number(atom);
            if (r_group_num) {
                attachment_point_nums.push_back(*r_group_num);
                BOOST_TEST(atom->getAtomicNum() == 0);
                BOOST_TEST(atom->getDegree() == 1);
            }
            BOOST_TEST(
                !atom->hasProp(RDKit::common_properties::molAtomMapNumber));
        }
        std::ranges::sort(attachment_point_nums);
        const std::vector<unsigned int> expected_nums = {1, 2};
        BOOST_TEST(attachment_point_nums == expected_nums);
    }
}

BOOST_AUTO_TEST_CASE(test_get_attachment_points_for_res)
{
    const std::vector<std::pair<int, std::string>> expected = {{1, "N"},
                                                               {2, "O"}};
    BOOST_TEST(get_attachment_points_for_res(
                   "A", rdkit_extensions::ChainType::PEPTIDE) == expected);
}

BOOST_AUTO_TEST_CASE(test_get_attachment_points_for_unknown_peptide)
{
    const std::vector<std::pair<int, std::string>> expected = {{1, ""},
                                                               {2, ""}};
    BOOST_TEST(get_attachment_points_for_res(
                   "X", rdkit_extensions::ChainType::PEPTIDE) == expected);
}

BOOST_AUTO_TEST_CASE(test_get_attachment_points_for_unknown_nucleic_acid)
{
    using rdkit_extensions::ChainType;
    const std::vector<std::pair<int, std::string>> expected = {{1, ""}};
    BOOST_TEST(get_attachment_points_for_res("N", ChainType::RNA) == expected);
}

BOOST_AUTO_TEST_CASE(test_get_attachment_points_for_missing_monomer)
{
    using rdkit_extensions::ChainType;
    const std::vector<std::pair<int, std::string>> expected = {{1, ""},
                                                               {2, ""}};
    for (const auto chain_type :
         {ChainType::PEPTIDE, ChainType::RNA, ChainType::CHEM}) {
        BOOST_TEST(get_attachment_points_for_res("missingMonomer",
                                                 chain_type) == expected);
    }
    // The unknown-base special case must not apply to other polymer types.
    BOOST_TEST(get_attachment_points_for_res("N", ChainType::CHEM) == expected);
}

/**
 * Make sure that contains_two_monomer_linkages correctly detects two monomer
 * linkages in the same bond when there's a disulfide bond between neighboring
 * cysteines.
 */
BOOST_AUTO_TEST_CASE(test_contains_two_monomer_linkages)
{
    // two neighboring cysteines, but no disulfide
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{C.C}$$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(0)));

    // two neighboring cysteines with a disulfide
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(contains_two_monomer_linkages(mol->getBondWithIdx(0)));

    // a disulfide, but between two non-neighboring cysteines
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0");
    BOOST_TEST(mol->getNumBonds() == 3);
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(0)));
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(1)));
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(2)));
}

/**
 * An isolated custom CHEM monomer must expose its first attachment point
 * without trying to find the highest numbered bound attachment point.
 */
BOOST_AUTO_TEST_CASE(test_isolated_chem_attachment_points)
{
    auto mol = rdkit_extensions::to_rdkit("CHEM1{[[*:1]C]}$$$$V2.0");
    prepare_mol(*mol);
    BOOST_REQUIRE(mol->getNumAtoms() == 1);
    BOOST_REQUIRE(mol->getNumBonds() == 0);

    auto [bound_aps, unbound_aps] =
        get_attachment_points_for_monomer(mol->getAtomWithIdx(0));
    const std::vector<UnboundAttachmentPoint> expected = {
        {"R1", "R1", 1, Direction::W}};
    BOOST_TEST(bound_aps.empty());
    BOOST_TEST(unbound_aps == expected);
}

/**
 * Make sure that get_bound_attachment_point_names_and_atoms() and
 * get_available_attachment_point_names() return the expected attachment point
 * names for a variety of molecules
 */
BOOST_AUTO_TEST_CASE(test_get_attachment_points)
{
    std::vector<BoundAttachmentPoint> bound_aps, exp_bound;
    std::vector<UnboundAttachmentPoint> unbound_aps, exp_available;
    RDKit::Atom *atom0, *atom1, *atom2;
    auto set_smiles_monomer = [](RDKit::Atom* atom, const std::string& smiles) {
        atom->setProp(SMILES_MONOMER, true);
        atom->setProp(ATOM_LABEL, smiles);
    };

    // a lone alanine has no bound attachment points
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"R1", "N", 1, Direction::W},
                         {"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // two alanines next to each other
    mol = rdkit_extensions::to_rdkit("PEPTIDE1{A.A}$$$$V2.0");
    prepare_mol(*mol);
    {
        const auto* atom0 = mol->getAtomWithIdx(0);
        const auto* atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // a side-chain interaction between two non-adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0");
    prepare_mol(*mol);
    {
        // put the three residues in a horizontal line
        auto& conf = mol->getConformer();
        conf.setAtomPos(0, {-BOND_LENGTH, 0.0, 0.0});
        conf.setAtomPos(1, {0.0, 0.0, 0.0});
        conf.setAtomPos(2, {BOND_LENGTH, 0.0, 0.0});

        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);
        atom2 = mol->getAtomWithIdx(2);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R3", "", 3, atom2, false, Direction::N}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R2", "C", 2, atom2, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom2);
        exp_bound = {{"R1", "N", 1, atom1, false, Direction::W},
                     {"R3", "", 3, atom0, false, Direction::N}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // a side-chain interaction between two adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R3", "", 3, atom1, true, Direction::N}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R3", "", 3, atom0, true, Direction::N}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // CHEM monomers
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[C([CH3:1])[CH3:2]]}|"
        "CHEM2{[C([CH3:1])([CH3:2])([CH3:3])[CH3:4]]}"
        "$CHEM1,CHEM2,1:R1-1:R3$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);
        set_smiles_monomer(atom0, "C([CH3:1])[CH3:2]");
        set_smiles_monomer(atom1, "C([CH3:1])([CH3:2])([CH3:3])[CH3:4]");

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", "R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", "R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R3", "R3", 3, atom0, false, Direction::N}};
        exp_available = {{"R1", "R1", 1, Direction::W},
                         {"R2", "R2", 2, Direction::E},
                         {"R4", "R4", 4, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // a bound attachment point absent from the CHEM definition doesn't invent
    // additional unbound attachment points
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[C([CH3:1])[CH3:2]]}|"
        "CHEM2{[C([CH3:1])([CH3:2])([CH3:3])[CH3:4]]}"
        "$CHEM1,CHEM2,1:R1-1:R11$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);
        set_smiles_monomer(atom0, "C([CH3:1])[CH3:2]");
        set_smiles_monomer(atom1, "C([CH3:1])([CH3:2])([CH3:3])[CH3:4]");

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", "R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", "R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R11", "R11", 11, atom0, false, Direction::N}};
        exp_available = {{"R1", "R1", 1, Direction::W},
                         {"R2", "R2", 2, Direction::E},
                         {"R3", "R3", 3, Direction::S},
                         {"R4", "R4", 4, Direction::NW}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // peptide R3 names reflect whether the attachment-point atom is sulfur
    for (const auto& [side_chain, expected_name] :
         std::array<std::pair<std::string, std::string>, 2>{
             {{"C[*:3]", "X"}, {"S[*:3]", "S"}}}) {
        const auto smiles =
            fmt::format("N([*:1])[C@@H]({})C(=O)[*:2]", side_chain);
        const auto helm = fmt::format("PEPTIDE1{{[{}]}}$$$$V2.0", smiles);
        mol = rdkit_extensions::to_rdkit(helm);
        prepare_mol(*mol);

        atom0 = mol->getAtomWithIdx(0);
        set_smiles_monomer(atom0, smiles);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"R1", "N", 1, Direction::W},
                         {"R2", "C", 2, Direction::E},
                         {"R3", expected_name, 3, Direction::N}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // NA_PHOSPHATE monomers name their unbound attachment points after the
    // attachment point of the bound sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P}$$$$");
    prepare_mol(*mol);
    {
        auto start_phos = mol->getAtomWithIdx(0);
        auto start_sugar = mol->getAtomWithIdx(1);
        auto middle_phos = mol->getAtomWithIdx(3);
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phosphate = mol->getAtomWithIdx(6);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(start_phos);
        exp_bound = {{"R2", "", 2, start_sugar, false, Direction::E}};
        exp_available = {{"R1", "5'", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(middle_phos);
        exp_bound = {{"R1", "", 1, start_sugar, false, Direction::W},
                     {"R2", "", 2, term_sugar, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phosphate);
        exp_bound = {{"R1", "", 1, term_sugar, false, Direction::W}};
        exp_available = {{"R2", "3'", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // NA_PHOSPHATE monomers still take their name from the bound sugar even if
    // they're at the end of a chain of phosphates
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P.P.P}$$$$");
    prepare_mol(*mol);
    {
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phos_chain_1 = mol->getAtomWithIdx(6);
        auto term_phos_chain_2 = mol->getAtomWithIdx(7);
        auto term_phos_chain_3 = mol->getAtomWithIdx(8);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_1);
        exp_bound = {{"R1", "", 1, term_sugar, false, Direction::W},
                     {"R2", "", 2, term_phos_chain_2, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_2);
        exp_bound = {{"R1", "", 1, term_phos_chain_1, false, Direction::W},
                     {"R2", "", 2, term_phos_chain_3, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_3);
        exp_bound = {{"R1", "", 1, term_phos_chain_2, false, Direction::W}};
        exp_available = {{"R2", "3'", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // attachment points on lone phosphates are unnamed since there's no bound
    // sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P}$$$$");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"R1", "", 1, Direction::W},
                         {"R2", "", 2, Direction::E}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // an amino acid with unrecognized attachment points
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{A.A}$PEPTIDE1,PEPTIDE1,1:R4-2:R4$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R2", "C", 2, atom1, false, Direction::E},
                     {"R4", "R4", 4, atom1, true, Direction::E}};
        exp_available = {{"R1", "N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R1", "N", 1, atom0, false, Direction::W},
                     {"R4", "R4", 4, atom0, true, Direction::W}};
        exp_available = {{"R2", "C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // single stranded RNA
    mol = rdkit_extensions::to_rdkit("RNA1{R(A)P.R(C)P.R(G)P}$$$$V2.0");
    prepare_mol(*mol);
    {
        auto sugar = mol->getAtomWithIdx(0);
        auto base = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(base);
        exp_bound = {{"R1", "N1/9", 1, sugar, false, Direction::N}};
        exp_available = {{"pair", "H-bond", ATTACHMENT_POINT_WITH_CUSTOM_NAME,
                          Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // double stranded RNA
    mol = rdkit_extensions::to_rdkit(
        "RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}"
        "$RNA1,RNA2,2:pair-9:pair"
        "|RNA1,RNA2,5:pair-6:pair"
        "|RNA1,RNA2,8:pair-3:pair$$$V2.0");
    prepare_mol(*mol);
    {
        auto sugar = mol->getAtomWithIdx(0);
        auto base = mol->getAtomWithIdx(1);
        auto paired_base = mol->getAtomWithIdx(17);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(base);
        exp_bound = {{"R1", "N1/9", 1, sugar, false, Direction::N},
                     {"pair", "pair", ATTACHMENT_POINT_WITH_CUSTOM_NAME,
                      paired_base, false, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());
    }
}

/**
 * Make sure that get_attachment_points returns the correct number for parseable
 * chain name, and returns -1 when the chain name can't be parsed.
 */
BOOST_AUTO_TEST_CASE(test_get_chain_num)
{
    using rdkit_extensions::ChainType;
    BOOST_TEST(get_chain_num("PEPTIDE1", ChainType::PEPTIDE) == 1);
    BOOST_TEST(get_chain_num("PEPTIDE3", ChainType::PEPTIDE) == 3);
    BOOST_TEST(get_chain_num("PEPTIDE", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("ABCDEFG", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("ABCDEFG2", ChainType::PEPTIDE) == -1);
    BOOST_TEST(get_chain_num("PEPTIDE1", ChainType::RNA) == -1);
    BOOST_TEST(get_chain_num("RNA2", ChainType::RNA) == 2);
}

/**
 * Make sure that get_first_available_chain_name returns the expected chain name
 * of the appropriate chain type.
 */
BOOST_AUTO_TEST_CASE(test_get_first_available_chain_name)
{
    using rdkit_extensions::ChainType;

    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}|PEPTIDE2{A}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE3");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA1");

    mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}|PEPTIDE3{A}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE2");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA1");

    mol = rdkit_extensions::to_rdkit(
        "RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}$$$$V2.0");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::RNA) == "RNA3");
    BOOST_TEST(get_first_available_chain_name(*mol, ChainType::PEPTIDE) ==
               "PEPTIDE1");
}

BOOST_AUTO_TEST_CASE(test_peptide_has_ap3)
{
    BOOST_TEST(peptide_has_ap3("A", false) == false);
    BOOST_TEST(peptide_has_ap3("K", false) == true);
    BOOST_TEST(peptide_has_ap3("C", false) == true);
    BOOST_TEST(peptide_has_ap3("CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2]", true) ==
               false);
    BOOST_TEST(peptide_has_ap3("O=C([C@H](CCCCN[H:3])N[H:1])[OH:2]", true) ==
               true);
    BOOST_TEST(peptide_has_ap3("O=C([C@H](CS[H:3])N[H:1])[OH:2]", true) ==
               true);
}

} // namespace sketcher
} // namespace schrodinger
