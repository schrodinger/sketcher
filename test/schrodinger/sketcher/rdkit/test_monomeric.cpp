
#define BOOST_TEST_MODULE monomeric

#include <unordered_set>

#include <rdkit/GraphMol/RWMol.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "schrodinger/rdkit_extensions/convert.h"
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
    std::vector<std::pair<std::string, const RDKit::Atom*>> exp_bound;
    std::vector<std::string> exp_available;
    const RDKit::Atom *atom0, *atom1, *atom2;

    // a lone alanine has no bound attachment points
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    {
        atom0 = mol->getAtomWithIdx(0);
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0).empty());
        exp_available = {"N", "C", "X"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);
    }

    // two alanines next to each other
    mol = rdkit_extensions::to_rdkit("PEPTIDE1{A.A}$$$$V2.0");
    {
        const auto* atom0 = mol->getAtomWithIdx(0);
        const auto* atom1 = mol->getAtomWithIdx(1);

        exp_bound = {{"C", atom1}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0) ==
                   exp_bound);
        exp_available = {"N", "X"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);

        exp_bound = {{"N", atom0}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom1) ==
                   exp_bound);
        exp_available = {"C", "X"};
        BOOST_TEST(get_available_attachment_point_names(atom1) ==
                   exp_available);
    }

    // a side-chain interaction between two non-adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0");
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);
        atom2 = mol->getAtomWithIdx(2);

        exp_bound = {{"C", atom1}, {"X", atom2}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0) ==
                   exp_bound);
        exp_available = {"N"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);

        exp_bound = {{"N", atom0}, {"C", atom2}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom1) ==
                   exp_bound);
        exp_available = {"X"};
        BOOST_TEST(get_available_attachment_point_names(atom1) ==
                   exp_available);

        exp_bound = {{"N", atom1}, {"X", atom0}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom2) ==
                   exp_bound);
        exp_available = {"C"};
        BOOST_TEST(get_available_attachment_point_names(atom2) ==
                   exp_available);
    }

    // a side-chain interaction between two adjacent residues
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        exp_bound = {{"C", atom1}, {"X", atom1}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0) ==
                   exp_bound);
        exp_available = {"N"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);

        exp_bound = {{"N", atom0}, {"X", atom0}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom1) ==
                   exp_bound);
        exp_available = {"C"};
        BOOST_TEST(get_available_attachment_point_names(atom1) ==
                   exp_available);
    }

    // CHEM monomers
    mol = rdkit_extensions::to_rdkit(
        "CHEM1{[MONO1]}|CHEM2{[MONO2]}$CHEM1,CHEM2,1:R1-1:R3$$$V2.0");
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        exp_bound = {{"R1", atom1}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0) ==
                   exp_bound);
        exp_available = {"R2"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);

        exp_bound = {{"R3", atom0}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom1) ==
                   exp_bound);
        exp_available = {"R1", "R2", "R4"};
        BOOST_TEST(get_available_attachment_point_names(atom1) ==
                   exp_available);
    }

    // NA_PHOSPHATE monomers name their unbound attachment points after the
    // attachment point of the bound sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P}$$$$");
    {
        auto start_phos = mol->getAtomWithIdx(0);
        auto start_sugar = mol->getAtomWithIdx(1);
        auto middle_phos = mol->getAtomWithIdx(3);
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phosphate = mol->getAtomWithIdx(6);

        exp_bound = {{"", start_sugar}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(start_phos) ==
                   exp_bound);
        exp_available = {"3'"};
        BOOST_TEST(get_available_attachment_point_names(start_phos) ==
                   exp_available);

        exp_bound = {{"", start_sugar}, {"", term_sugar}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(middle_phos) ==
                   exp_bound);
        BOOST_TEST(get_available_attachment_point_names(middle_phos).empty());
        
        
        exp_bound = {{"", term_sugar}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(term_phosphate) ==
                   exp_bound);
        exp_available = {"5'"};
        BOOST_TEST(get_available_attachment_point_names(term_phosphate) ==
                   exp_available);
    }

    // NA_PHOSPHATE monomers still take their name from the bound sugar even if
    // they're at the end of a chain of phosphates
    mol = rdkit_extensions::to_rdkit("RNA1{P.R(U)P.R(T)P.P.P}$$$$");
    {
        auto term_sugar = mol->getAtomWithIdx(4);
        auto term_phos_chain_1 = mol->getAtomWithIdx(6);
        auto term_phos_chain_2 = mol->getAtomWithIdx(7);
        auto term_phos_chain_3 = mol->getAtomWithIdx(8);

        exp_bound = {{"", term_sugar}, {"", term_phos_chain_2}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(
                       term_phos_chain_1) == exp_bound);
        BOOST_TEST(
            get_available_attachment_point_names(term_phos_chain_1).empty());

        exp_bound = {{"", term_phos_chain_1}, {"", term_phos_chain_3}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(
                       term_phos_chain_2) == exp_bound);
        BOOST_TEST(
            get_available_attachment_point_names(term_phos_chain_2).empty());

        exp_bound = {{"", term_phos_chain_2}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(
                       term_phos_chain_3) == exp_bound);
        exp_available = {"5'"};
        BOOST_TEST(get_available_attachment_point_names(term_phos_chain_3) ==
                   exp_available);
    }

    // attachment points on lone phosphates are unnamed since there's no bound
    // sugar
    mol = rdkit_extensions::to_rdkit("RNA1{P}$$$$");
    {
        atom0 = mol->getAtomWithIdx(0);
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0).empty());
        exp_available = {"", ""};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);
    }

    // an amino acid with an unrecognized attachment point (R4)
    mol = rdkit_extensions::to_rdkit(
        "PEPTIDE1{A.A}$PEPTIDE1,PEPTIDE1,1:R3-2:R4$$$V2.0");
    {
        atom0 = mol->getAtomWithIdx(0);
        atom1 = mol->getAtomWithIdx(1);

        exp_bound = {{"C", atom1}, {"X", atom1}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom0) ==
                   exp_bound);
        exp_available = {"N"};
        BOOST_TEST(get_available_attachment_point_names(atom0) ==
                   exp_available);

        exp_bound = {{"N", atom0}, {"R4", atom0}};
        BOOST_TEST(get_bound_attachment_point_names_and_atoms(atom1) ==
                   exp_bound);
        exp_available = {"C", "X"};
        BOOST_TEST(get_available_attachment_point_names(atom1) ==
                   exp_available);
    }
}

} // namespace sketcher
} // namespace schrodinger
