#define BOOST_TEST_MODULE test_draw_nucleotide_scene_tool_integration_tests

#include <algorithm>
#include <vector>

#include <QPointF>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/helm.h"
#include "./test_monomer_integration_tests_common.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{
// SKETCH-2765 helper: gather every sugar->base (branch) bond length in the
// current mol so the test can assert they stay uniform as the chain grows.
std::vector<double> collect_branch_bond_lengths(MolModel* mol_model)
{
    auto mol = mol_model->getMol();
    BOOST_REQUIRE(mol);
    const auto& conf = mol->getConformer();
    std::vector<double> lengths;
    for (const auto* bond : mol->bonds()) {
        const auto* a = bond->getBeginAtom();
        const auto* b = bond->getEndAtom();
        bool a_branch =
            a->hasProp(BRANCH_MONOMER) && a->getProp<bool>(BRANCH_MONOMER);
        bool b_branch =
            b->hasProp(BRANCH_MONOMER) && b->getProp<bool>(BRANCH_MONOMER);
        if (!a_branch && !b_branch) {
            continue;
        }
        const auto& pa = conf.getAtomPos(a->getIdx());
        const auto& pb = conf.getAtomPos(b->getIdx());
        lengths.push_back((pa - pb).length());
    }
    return lengths;
}

// Canonical sugar->base distance produced by coordgen
// (rdkit_extensions::MONOMER_BOND_LENGTH = SIDE_TO_SIDE_DISTANCE 0.7 +
// MONOMER_MINIMUM_SIZE 0.8).
constexpr double EXPECTED_BRANCH_BOND_LENGTH = 1.5;

void check_branch_bonds_stable(MolModel* mol_model, size_t expected_count,
                               int add_number)
{
    auto lengths = collect_branch_bond_lengths(mol_model);
    BOOST_REQUIRE_EQUAL(lengths.size(), expected_count);
    auto [min_it, max_it] = std::minmax_element(lengths.begin(), lengths.end());
    BOOST_TEST(*max_it - *min_it < 1e-6,
               "branch bond lengths drifted after add #"
                   << add_number << ": min=" << *min_it << " max=" << *max_it);
    // Pin to the canonical value so a future bug that uniformly shrinks or
    // stretches every bond can't slip past the uniformity check.
    BOOST_TEST(std::abs(*min_it - EXPECTED_BRANCH_BOND_LENGTH) < 1e-6,
               "branch bond length is " << *min_it << " after add #"
                                        << add_number << ", expected "
                                        << EXPECTED_BRANCH_BOND_LENGTH);
}
} // namespace

// SKETCH-2765: adding R(A)P units used to drift the sugar->base bond longer
// on every add — resize_monomers treated first-time sizing of a new monomer
// as a resize event and pushed every other monomer outward.
BOOST_AUTO_TEST_CASE(test_sketch_2765_branch_bond_growth_separate_polymers)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);

    for (int i = 0; i < 7; ++i) {
        QPointF pos(i * 200.0, 0.0);
        fix.mouseClick(pos);
        check_branch_bonds_stable(fix.m_mol_model, static_cast<size_t>(i + 1),
                                  i + 1);
    }
}

// SKETCH-2765: same regression, but exercised against a single connected chain
// (the actual user-reported repro from the screenshot).
BOOST_AUTO_TEST_CASE(test_sketch_2765_branch_bond_growth_single_chain)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});

    for (int i = 1; i < 7; ++i) {
        // The RNA nucleotide tool emits monomers in (R, A, P) triples in that
        // canonical HELM order, so the trailing phosphate of the i-th
        // nucleotide sits at atom index 3*i - 1 (2, 5, 8, ...).
        auto phos_pos = fix.getMonomerPos(3 * i - 1);
        fix.mouseMove(phos_pos);
        fix.mouseClick(phos_pos);
        check_branch_bonds_stable(fix.m_mol_model, static_cast<size_t>(i + 1),
                                  i + 1);
    }
}

BOOST_AUTO_TEST_CASE(test_click_empty_space_adds_nucleotide)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{R(A)P}$$$$V2.0");

    fix.setRNANucleotideTool(StdNucleobase::U_OR_T);
    fix.mouseClick({200, 0});
    fix.verifyHELM("RNA1{R(A)P}|RNA2{R(U)P}$$$$V2.0");

    // create a DNA nucleotide
    fix.setDNANucleotideTool(StdNucleobase::U_OR_T);
    fix.mouseClick({400, 0});
    fix.verifyHELM("RNA1{R(A)P}|RNA2{R(U)P}|RNA3{[dR](T)P}$$$$V2.0");

    // create a custom nucleotide
    fix.setCustomNucleotideTool("Tho", "I", "PS");
    fix.mouseClick({600, 0});
    fix.verifyHELM(
        "RNA1{R(A)P}|RNA2{R(U)P}|RNA3{[dR](T)P}|RNA4{[Tho](I)[PS]}$$$$V2.0");
}

BOOST_AUTO_TEST_CASE(test_click_existing_monomers)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{R(A)P}$$$$V2.0");

    // click on the phosphate of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::C);
    auto phos_pos = fix.getMonomerPos(2);
    fix.mouseMove(phos_pos);
    fix.mouseClick(phos_pos);
    fix.verifyHELM("RNA1{R(A)P.R(C)P}$$$$V2.0");

    // click on the sugar of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::G);
    auto sugar_pos = fix.getMonomerPos(0);
    fix.mouseMove(sugar_pos);
    fix.mouseClick(sugar_pos);
    fix.verifyHELM("RNA1{R(G)P.R(A)P.R(C)P}$$$$V2.0");

    // click on the base of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::U_OR_T);
    auto base_pos = fix.getMonomerPos(1);
    fix.mouseMove(base_pos);
    fix.mouseClick(base_pos);
    fix.verifyHELM(
        "RNA1{R(G)P.R(A)P.R(C)P}|RNA2{R(U)P}$RNA1,RNA2,5:pair-2:pair$$$V2.0");

    // clicking on the same spots again shouldn't have any effect since those
    // monomers no longer have attachment points available

    // recalculate the coordinates first in case anything has been adjusted to
    // fit the labels
    sugar_pos = fix.getMonomerPos(0);
    base_pos = fix.getMonomerPos(1);
    phos_pos = fix.getMonomerPos(2);

    fix.mouseMove(sugar_pos);
    fix.mouseClick(sugar_pos);
    fix.verifyHELM(
        "RNA1{R(G)P.R(A)P.R(C)P}|RNA2{R(U)P}$RNA1,RNA2,5:pair-2:pair$$$V2.0");

    fix.mouseMove(phos_pos);
    fix.mouseClick(phos_pos);
    fix.verifyHELM(
        "RNA1{R(G)P.R(A)P.R(C)P}|RNA2{R(U)P}$RNA1,RNA2,5:pair-2:pair$$$V2.0");
    fix.mouseMove(base_pos);
    fix.mouseClick(base_pos);
    fix.verifyHELM(
        "RNA1{R(G)P.R(A)P.R(C)P}|RNA2{R(U)P}$RNA1,RNA2,5:pair-2:pair$$$V2.0");
}

BOOST_AUTO_TEST_CASE(test_click_existing_monomer_attachment_points)
{
    MonomerToolTestFixture fix;
    fix.setNucleicAcidTool(NucleicAcidTool::R);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{R}$$$$V2.0");

    fix.setRNANucleotideTool(StdNucleobase::A);
    // hover over the ribose to make the attachment points appear
    fix.mouseMove({0, 0});

    // clicking on the 3' or 1' attachment points should have no effect
    auto ap_3p_pos = fix.getAttachmentPointPos(0, "3'");
    fix.mouseMove(ap_3p_pos);
    fix.mouseClick(ap_3p_pos);
    fix.verifyHELM("RNA1{R}$$$$V2.0");

    auto ap_1p_pos = fix.getAttachmentPointPos(0, "1'");
    fix.mouseMove(ap_1p_pos);
    fix.mouseClick(ap_1p_pos);
    fix.verifyHELM("RNA1{R}$$$$V2.0");

    // clicking on the 5' attachment point should add a nucleotide
    auto ap_5p_pos = fix.getAttachmentPointPos(0, "5'");
    fix.mouseMove(ap_5p_pos);
    fix.mouseClick(ap_5p_pos);
    fix.verifyHELM("RNA1{R(A)P.R}$$$$V2.0");
}

// Regression test: dragging from an existing nucleotide bead into empty space
// and then moving the mouse used to crash because onLeftButtonRelease deleted
// m_hint_fragment_item without nulling it, leaving a dangling pointer that the
// next mouse-move's clearHintFragmentItem() would double-free.
BOOST_AUTO_TEST_CASE(test_drag_from_bead_to_empty_then_move_does_not_crash)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{R(A)P}$$$$V2.0");

    // Re-arm the nucleotide tool so the scene tool is the fragment tool again
    fix.setRNANucleotideTool(StdNucleobase::C);

    // Hover over the sugar bead so the fragment hint is created
    auto sugar_pos = fix.getMonomerPos(0);
    fix.mouseMove(sugar_pos);

    // Press on the bead, drag to empty space, release
    QPointF empty_pos{500, 500};
    fix.mousePress(sugar_pos);
    fix.mouseMove(empty_pos, Qt::LeftButton);
    fix.mouseRelease(empty_pos);

    // Moving the mouse after the release used to dereference the freed
    // m_hint_fragment_item and crash. The structure must be unchanged.
    fix.mouseMove(empty_pos + QPointF(1, 1));
    fix.verifyHELM("RNA1{R(A)P}$$$$V2.0");
}

} // namespace sketcher
} // namespace schrodinger
