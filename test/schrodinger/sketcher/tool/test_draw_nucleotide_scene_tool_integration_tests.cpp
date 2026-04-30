#define BOOST_TEST_MODULE test_draw_nucleotide_scene_tool_integration_tests

#include <QPointF>
#include <boost/test/unit_test.hpp>

#include "./test_monomer_integration_tests_common.h"

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_click_empty_space_adds_nucleotide)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{P.R(A)}$$$$V2.0");

    fix.setRNANucleotideTool(StdNucleobase::U_OR_T);
    fix.mouseClick({200, 0});
    fix.verifyHELM("RNA1{P.R(A)}|RNA2{P.R(U)}$$$$V2.0");

    // create a DNA nucleotide
    fix.setDNANucleotideTool(StdNucleobase::U_OR_T);
    fix.mouseClick({400, 0});
    fix.verifyHELM("RNA1{P.R(A)}|RNA2{P.R(U)}|RNA3{P.[dR](T)}$$$$V2.0");

    // create a custom nucleotide
    fix.setCustomNucleotideTool("Tho", "I", "PS");
    fix.mouseClick({600, 0});
    fix.verifyHELM("RNA1{P.R(A)}|RNA2{P.R(U)}|RNA3{P.[dR](T)}|RNA4{[PS].[Tho]("
                   "I)}$$$$V2.0");
}

BOOST_AUTO_TEST_CASE(test_click_existing_monomers)
{
    MonomerToolTestFixture fix;
    fix.setRNANucleotideTool(StdNucleobase::A);
    fix.mouseClick({0, 0});
    fix.verifyHELM("RNA1{P.R(A)}$$$$V2.0");

    // click on the sugar of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::C);
    auto sugar_pos = fix.getMonomerPos(1);
    fix.mouseMove(sugar_pos);
    fix.mouseClick(sugar_pos);
    fix.verifyHELM("RNA1{P.R(A)P.R(C)}$$$$V2.0");

    // click on the phosphate of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::G);
    auto phos_pos = fix.getMonomerPos(0);
    fix.mouseMove(phos_pos);
    fix.mouseClick(phos_pos);
    fix.verifyHELM("RNA1{P.R(G)P.R(A)P.R(C)}$$$$V2.0");

    // click on the base of the first nucleotide
    fix.setRNANucleotideTool(StdNucleobase::U_OR_T);
    auto base_pos = fix.getMonomerPos(2);
    fix.mouseMove(base_pos);
    fix.mouseClick(base_pos);
    fix.verifyHELM(
        "RNA1{P.R(G)P.R(A)P.R(C)}|RNA2{P.R(U)}$RNA1,RNA2,6:pair-3:pair$$$V2.0");

    // clicking on the same spots again shouldn't have any effect since those
    // monomers no longer have attachment points available

    // recalculate the coordinates first in case anything has been adjusted to
    // fit the labels
    sugar_pos = fix.getMonomerPos(1);
    phos_pos = fix.getMonomerPos(0);
    base_pos = fix.getMonomerPos(2);

    fix.mouseMove(sugar_pos);
    fix.mouseClick(sugar_pos);
    fix.verifyHELM(
        "RNA1{P.R(G)P.R(A)P.R(C)}|RNA2{P.R(U)}$RNA1,RNA2,6:pair-3:pair$$$V2.0");

    fix.mouseMove(phos_pos);
    fix.mouseClick(phos_pos);
    fix.verifyHELM(
        "RNA1{P.R(G)P.R(A)P.R(C)}|RNA2{P.R(U)}$RNA1,RNA2,6:pair-3:pair$$$V2.0");

    fix.mouseMove(base_pos);
    fix.mouseClick(base_pos);
    fix.verifyHELM(
        "RNA1{P.R(G)P.R(A)P.R(C)}|RNA2{P.R(U)}$RNA1,RNA2,6:pair-3:pair$$$V2.0");
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
    fix.verifyHELM("RNA1{R.P.R(A)}$$$$V2.0");
}

} // namespace sketcher
} // namespace schrodinger
