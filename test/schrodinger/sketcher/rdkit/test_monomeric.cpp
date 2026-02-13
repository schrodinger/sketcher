
#define BOOST_TEST_MODULE monomeric

#include <unordered_set>

#include <rdkit/GraphMol/RWMol.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

using namespace boost::unit_test;

namespace schrodinger
{
namespace sketcher
{

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
 * Make sure that get_bound_attachment_point_names_and_atoms() and
 * get_available_attachment_point_names() return the expected attachment point
 * names for a variety of molecules
 */
BOOST_AUTO_TEST_CASE(test_get_attachment_points)
{
    std::vector<BoundAttachmentPoint> bound_aps, exp_bound;
    std::vector<UnboundAttachmentPoint> unbound_aps, exp_available;
    const RDKit::Atom *atom0, *atom1, *atom2;

    // a lone alanine has no bound attachment points
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_available = {{"N", 1, Direction::W},
                         {"C", 2, Direction::E},
                         {"X", 3, Direction::N}};
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
        exp_bound = {{"C", 2, atom1, false, Direction::E}};
        exp_available = {{"N", 1, Direction::W}, {"X", 3, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"N", 1, atom0, false, Direction::W}};
        exp_available = {{"C", 2, Direction::E}, {"X", 3, Direction::N}};
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
        exp_bound = {{"C", 2, atom1, false, Direction::E},
                     {"X", 3, atom2, false, Direction::N}};
        exp_available = {{"N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"N", 1, atom0, false, Direction::W},
                     {"C", 2, atom2, false, Direction::E}};
        exp_available = {{"X", 3, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom2);
        exp_bound = {{"N", 1, atom1, false, Direction::W},
                     {"X", 3, atom0, false, Direction::N}};
        exp_available = {{"C", 2, Direction::E}};
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
        exp_bound = {{"C", 2, atom1, false, Direction::E},
                     {"X", 3, atom1, true, Direction::N}};
        exp_available = {{"N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"N", 1, atom0, false, Direction::W},
                     {"X", 3, atom0, true, Direction::N}};
        exp_available = {{"C", 2, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // CHEM monomers
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[MONO1]}|CHEM2{[MONO2]}$CHEM1,CHEM2,1:R1-1:R3$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R3", 3, atom0, false, Direction::N}};
        exp_available = {{"R1", 1, Direction::W},
                         {"R2", 2, Direction::E},
                         {"R4", 4, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);
    }

    // CHEM monomers with too many attachment points - some will be omitted
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[MONO1]}|CHEM2{[MONO2]}$CHEM1,CHEM2,1:R1-1:R11$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"R1", 1, atom1, false, Direction::S}};
        exp_available = {{"R2", 2, Direction::N}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"R11", 11, atom0, false, Direction::N}};
        exp_available = {{"R1", 1, Direction::W},  {"R2", 2, Direction::E},
                         {"R3", 3, Direction::S},  {"R4", 4, Direction::NW},
                         {"R5", 5, Direction::NE}, {"R6", 6, Direction::SE},
                         {"R7", 7, Direction::SW}};
        BOOST_TEST(bound_aps == exp_bound);
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
        exp_bound = {{"", 2, start_sugar, false, Direction::E}};
        exp_available = {{"3'", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(middle_phos);
        exp_bound = {{"", 1, start_sugar, false, Direction::W},
                     {"", 2, term_sugar, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phosphate);
        exp_bound = {{"", 1, term_sugar, false, Direction::W}};
        exp_available = {{"5'", 2, Direction::E}};
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
        exp_bound = {{"", 1, term_sugar, false, Direction::W},
                     {"", 2, term_phos_chain_2, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_2);
        exp_bound = {{"", 1, term_phos_chain_1, false, Direction::W},
                     {"", 2, term_phos_chain_3, false, Direction::E}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(term_phos_chain_3);
        exp_bound = {{"", 1, term_phos_chain_2, false, Direction::W}};
        exp_available = {{"5'", 2, Direction::E}};
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
        exp_available = {{"", 1, Direction::W}, {"", 2, Direction::E}};
        BOOST_TEST(bound_aps.empty());
        BOOST_TEST(unbound_aps == exp_available);
    }

    // an amino acid with an unrecognized attachment point (R4)
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{A.A}$PEPTIDE1,PEPTIDE1,1:R3-2:R4$$$V2.0");
    prepare_mol(*mol);
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom0);
        exp_bound = {{"C", 2, atom1, false, Direction::E},
                     {"X", 3, atom1, true, Direction::E}};
        exp_available = {{"N", 1, Direction::W}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps == exp_available);

        std::tie(bound_aps, unbound_aps) =
            get_attachment_points_for_monomer(atom1);
        exp_bound = {{"N", 1, atom0, false, Direction::W},
                     {"R4", 4, atom0, true, Direction::W}};
        exp_available = {{"C", 2, Direction::E}, {"X", 3, Direction::N}};
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
        exp_bound = {{"N1/9", 1, sugar, false, Direction::N}};
        exp_available = {
            {"pair", ATTACHMENT_POINT_WITH_CUSTOM_NAME, Direction::S}};
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
        exp_bound = {{"N1/9", 1, sugar, false, Direction::N},
                     {"pair", ATTACHMENT_POINT_WITH_CUSTOM_NAME, paired_base,
                      false, Direction::S}};
        BOOST_TEST(bound_aps == exp_bound);
        BOOST_TEST(unbound_aps.empty());
    }
}

} // namespace sketcher
} // namespace schrodinger
