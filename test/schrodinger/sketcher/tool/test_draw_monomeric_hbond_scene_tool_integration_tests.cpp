#define BOOST_TEST_MODULE test_draw_monomeric_hbond_scene_tool_integration_tests

#include <QPointF>
#include <boost/test/unit_test.hpp>

#include "./test_monomer_integration_tests_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Confirm that clicking in empty space has no effect
 */
BOOST_AUTO_TEST_CASE(test_click_empty_space)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    fix.mouseClick({0, 0});
    fix.confirmIsEmpty();
}

/**
 * Confirm that clicking on an existing monomer has no effect
 */
BOOST_AUTO_TEST_CASE(test_click_existing_monomer)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto pos = fix.getMonomerPos(0);
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    fix.mouseClick(pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer to empty space has no
 * effect
 */
BOOST_AUTO_TEST_CASE(test_drag_monomer_to_empty)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);

    // Drag from the monomer to the right
    auto start_pos = fix.getMonomerPos(0);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer to an existing monomer
 * connects them via the default attachment points
 */
BOOST_AUTO_TEST_CASE(test_drag_monomer_to_monomer)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    fix.importMolText("PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{F}$$$$V2.0");
    auto pos1 = fix.getMonomerPos(0);
    auto pos2 = fix.getMonomerPos(1);
    auto pos3 = fix.getMonomerPos(2);
    fix.mouseDrag(pos1, pos2);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{F}$PEPTIDE1,PEPTIDE2,1:"
                   "pair-1:pair$$$V2.0");

    // drag between the same two monomers and make sure that we don't try to add
    // a second hydrogen bond
    fix.mouseDrag(pos1, pos2);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{F}$PEPTIDE1,PEPTIDE2,1:"
                   "pair-1:pair$$$V2.0");

    // create a hydrogen bond from chain 1 to chain 3
    fix.mouseDrag(pos1, pos3);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{F}$PEPTIDE1,PEPTIDE2,1:"
                   "pair-1:pair|PEPTIDE1,PEPTIDE3,1:pair-1:pair$$$V2.0");

    // create a hydrogen bond from chain 2 to chain 3
    fix.mouseDrag(pos2, pos3);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}|PEPTIDE3{F}$PEPTIDE1,PEPTIDE2,1:"
                   "pair-1:pair|PEPTIDE1,PEPTIDE3,1:pair-1:pair|PEPTIDE2,"
                   "PEPTIDE3,1:pair-1:pair$$$V2.0");
}

/**
 * Try dragging to and from a location where an attachment point would be with
 * the DrawMonomericConnectionSceneTool.  Since the DrawMonomericHBondSceneTool
 * doesn't display attachment points, these drags should have no effect.
 */
BOOST_AUTO_TEST_CASE(test_drag_to_and_from_ap_location)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    fix.importMolText("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    auto pos1 = fix.getMonomerPos(0);
    auto pos2 = fix.getMonomerPos(1);
    fix.mouseDrag(fix.getPosJustOutsideOfMonomer(pos1, Direction::E), pos2);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    fix.mouseDrag(fix.getPosJustOutsideOfMonomer(pos1, Direction::W), pos2);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    fix.mouseDrag(fix.getPosJustOutsideOfMonomer(pos1, Direction::N), pos2);
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    fix.mouseDrag(pos1, fix.getPosJustOutsideOfMonomer(pos2, Direction::E));
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    fix.mouseDrag(pos1, fix.getPosJustOutsideOfMonomer(pos2, Direction::W));
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    fix.mouseDrag(pos1, fix.getPosJustOutsideOfMonomer(pos2, Direction::N));
    fix.verifyHELM("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
}

/**
 * Confirm that dragging from empty space to empty space has no effect
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_empty)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_pos = QPointF(100, 100);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);
    fix.confirmIsEmpty();
}

/**
 * Confirm that dragging from empty space to an existing monomer has no effect
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_monomer)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{C}$$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_pos = QPointF(100, 100);
    auto end_pos = fix.getMonomerPos(0);
    fix.mouseDrag(start_pos, end_pos);
    fix.verifyHELM("PEPTIDE1{C}$$$$V2.0");
}

/**
 * Make sure that dragging between two monomers connected by a standard backbone
 * bond doesn't add a hydrogen bond between them, since monomer_mol won't mix
 * covalent and hydrogen bonds in a single connection.
 */
BOOST_AUTO_TEST_CASE(test_drag_already_backbone_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A.A}$$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{A.A}$$$$V2.0");
}

/**
 * Make sure that dragging between two monomers connected by a custom bond
 * doesn't add a hydrogen bond between them, since monomer_mol won't mix
 * covalent and hydrogen bonds in a single connection.
 */
BOOST_AUTO_TEST_CASE(test_drag_already_custom_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText(
        "PEPTIDE1{C}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,1:R2-1:R2$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM(
        "PEPTIDE1{C}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,1:R2-1:R2$$$V2.0");
}

/**
 * Make sure that dragging between two monomers connected by both backbone and
 * custom bonds doesn't add a hydrogen bond between them, since monomer_mol
 * won't mix covalent and hydrogen bonds in a single connection.
 */
BOOST_AUTO_TEST_CASE(test_drag_third_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
}

/**
 * Start a drag from one monomer toward another, then undo (which destroys the
 * second monomer) before releasing.
 */
BOOST_AUTO_TEST_CASE(test_undo_drag_end_monomer_no_crash)
{
    MonomerToolTestFixture fix;

    // place two monomers far enough apart that they don't overlap
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    fix.mouseClick(QPointF(100, 100));
    auto pos_a = fix.getMonomerPos(0);
    fix.mouseClick(pos_a + QPointF(200, 0));
    BOOST_REQUIRE(fix.m_mol_model->getMol()->getNumAtoms() == 2);
    auto pos_b = fix.getMonomerPos(1);

    // start a drag from A to B
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    fix.mouseMove(pos_a);
    fix.mousePress(pos_a);
    fix.mouseMove(pos_b, Qt::LeftButton);

    // undo the second click - destroys monomer B
    auto* undo_stack = fix.m_mol_model->findChild<QUndoStack*>();
    BOOST_REQUIRE(undo_stack != nullptr);
    undo_stack->undo();
    process_qt_events();
    BOOST_TEST(fix.m_mol_model->getMol()->getNumAtoms() == 1);

    // releasing the drag must not crash
    fix.mouseRelease(pos_b);
}

/**
 * Make sure that dragging from a nucleic acid base's "pair" attachment point to
 * a peptide monomer doesn't cause a crash.
 */
BOOST_AUTO_TEST_CASE(test_drag_from_na_base_to_peptide)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}|RNA1{A}$$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_monomer_pos = fix.getMonomerPos(1);
    auto end_monomer_pos = fix.getMonomerPos(0);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}|RNA1{A}$RNA1,PEPTIDE1,1:pair-1:pair$$$V2.0");
}

/**
 * Make sure that dragging to a nucleic acid base's "pair" attachment point from
 * a peptide monomer doesn't cause a crash (i.e. dragging in the opposite
 * direction relative to the last test).
 */
BOOST_AUTO_TEST_CASE(test_drag_to_na_base_from_peptide)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}|RNA1{A}$$$$V2.0");
    fix.setMonomericConnectionTool(MonomericConnectionTool::HBOND);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}|RNA1{A}$PEPTIDE1,RNA1,1:pair-1:pair$$$V2.0");
}

} // namespace sketcher
} // namespace schrodinger
