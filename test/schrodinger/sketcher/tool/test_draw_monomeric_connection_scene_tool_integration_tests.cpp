#define BOOST_TEST_MODULE \
    test_draw_monomeric_connection_scene_tool_integration_tests

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
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
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
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    fix.mouseClick(pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that clicking on an unbounnd attachment point of an existing monomer
 * has no effect
 */
BOOST_AUTO_TEST_CASE(test_click_attachment_point)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto monomer_pos = fix.getMonomerPos(0);
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    // hover over the monomer to trigger AP label creation
    fix.mouseMove(monomer_pos);

    // click on the N terminus attachment point
    fix.mouseMove(monomer_pos);
    auto n_ap_pos = fix.getAttachmentPointPos(0, "N");
    fix.mouseClick(n_ap_pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // click on the C terminus attachment point
    fix.mouseMove(monomer_pos);
    auto c_ap_pos = fix.getAttachmentPointPos(0, "C");
    fix.mouseClick(c_ap_pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // click on the side chain attachment point
    fix.mouseMove(monomer_pos);
    auto x_ap_pos = fix.getAttachmentPointPos(0, "X");
    fix.mouseClick(x_ap_pos);
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
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);

    // Drag from the monomer to the right
    auto start_pos = fix.getMonomerPos(0);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from the attachment point of an existing monomer
 * to empty space has no effect
 */
BOOST_AUTO_TEST_CASE(test_drag_ap_to_empty)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);

    // Add initial monomer
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");

    // hover over the monomer so that the attachment point graphics items are
    // created
    auto ala_pos = fix.getMonomerPos(0);
    fix.mouseMove(ala_pos);

    // Drag from N attachment point to empty space
    auto start_pos = fix.getAttachmentPointPos(0, "N");
    auto end_pos = start_pos + QPointF(-100, 0);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // hover over the first monomer so that its attachment point graphics items
    // are created again
    fix.mouseMove(end_pos);
    fix.mouseMove(ala_pos);

    // Drag from C attachment point to empty space
    start_pos = fix.getAttachmentPointPos(0, "C");
    end_pos = start_pos + QPointF(100, 100);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // hover over the first monomer so that its attachment point graphics items
    // are created again
    fix.mouseMove(end_pos);
    fix.mouseMove(ala_pos);

    // Drag from C attachment point to empty space
    start_pos = fix.getAttachmentPointPos(0, "X");
    end_pos = start_pos + QPointF(-50, 100);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer to an existing monomer
 * connects them via the default attachment points
 */
BOOST_AUTO_TEST_CASE(test_drag_monomer_to_monomer_connects_default_aps)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    fix.importMolText("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    auto pos1 = fix.getMonomerPos(0);
    auto pos2 = fix.getMonomerPos(1);
    fix.mouseDrag(pos1, pos2);
    fix.verifyHELM("PEPTIDE1{A.C}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from the attachment point of one existing monomer
 * to the attachment point of another existing monomer connects the monomer via
 * the specified attachment points.
 */
BOOST_AUTO_TEST_CASE(test_drag_ap_to_ap_connects_via_both_aps)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    fix.importMolText("PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0");
    auto ala_pos = fix.getMonomerPos(0);
    auto cys_pos = fix.getMonomerPos(1);
    fix.mouseMove(ala_pos);
    auto start_pos = fix.getAttachmentPointPos(0, "N");
    fix.mouseMove(start_pos);
    fix.mousePress(start_pos);
    // first, drag to the cysteine to make its attachment points appear
    fix.mouseMove(cys_pos, Qt::LeftButton);
    auto end_pos = fix.getAttachmentPointPos(1, "N");
    fix.mouseMove(end_pos, Qt::LeftButton);
    fix.mouseRelease(end_pos);
    fix.verifyHELM(
        "PEPTIDE1{A}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:R1-1:R1$$$V2.0");
}

/**
 * Confirm that dragging from empty space to empty space has no effect
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_empty)
{
    MonomerToolTestFixture fix;
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
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
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    auto start_pos = QPointF(100, 100);
    auto end_pos = fix.getMonomerPos(0);
    fix.mouseDrag(start_pos, end_pos);
    fix.verifyHELM("PEPTIDE1{C}$$$$V2.0");
}

/**
 * Make sure that dragging between two existing monomer won't create a second
 * standard backbone connection between them (since monomer_mol can't store more
 * than one standard backbone connection between the same pair of monomers)
 */
BOOST_AUTO_TEST_CASE(test_drag_second_backbone_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A.A}$$$$V2.0");
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(0, "N");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    auto end_ap_pos = fix.getAttachmentPointPos(1, "C");
    fix.mouseMove(end_ap_pos);
    fix.mouseRelease(end_ap_pos);
    fix.verifyHELM("PEPTIDE1{A.A}$$$$V2.0");
}

/**
 * Make sure that dragging between two existing monomer won't create a second
 * custom connection (i.e. not a standard backbone connection) between them
 * (since monomer_mol can't store more than one custom connection between the
 * same pair of monomers)
 */
BOOST_AUTO_TEST_CASE(test_drag_second_custom_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText(
        "PEPTIDE1{C}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,1:R2-1:R2$$$V2.0");
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(0, "X");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    auto end_ap_pos = fix.getAttachmentPointPos(1, "X");
    fix.mouseMove(end_ap_pos);
    fix.mouseRelease(end_ap_pos);
    fix.verifyHELM(
        "PEPTIDE1{C}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,1:R2-1:R2$$$V2.0");
}

/**
 * Make sure that dragging between two existing monomer won't create a third
 * connection between them (since monomer_mol can't store more than two
 * connections between the same pair of monomers)
 */
BOOST_AUTO_TEST_CASE(test_drag_third_connection)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(0, "N");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    auto end_ap_pos = fix.getAttachmentPointPos(1, "C");
    fix.mouseMove(end_ap_pos);
    fix.mouseRelease(end_ap_pos);
    fix.verifyHELM("PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
}

/**
 * Make sure that dragging between two existing monomer won't create a
 * connection between them if there's already a hydrogen bond (since monomer_mol
 * won't mix covalent and hydrogen bonds)
 */
BOOST_AUTO_TEST_CASE(test_drag_hydrogen_bond)
{
    MonomerToolTestFixture fix;
    fix.importMolText(
        "PEPTIDE1{A}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:pair-1:pair$$$V2.0");
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);

    // drag from monomer to monomer
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseDrag(start_monomer_pos, end_monomer_pos);
    fix.verifyHELM(
        "PEPTIDE1{A}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:pair-1:pair$$$V2.0");
}

/**
 * Place a monomer, hover next to it so the tool caches unbound attachment-
 * point items and a hint fragment, then undo the placement and move the
 * mouse. The monomer's child AP items get destroyed when the monomer is
 * removed; without onStructureUpdated() the tool's stale raw pointers cause
 * a double-free / use-after-free on the next mouse move.
 */
BOOST_AUTO_TEST_CASE(test_undo_monomer_after_hint_fragment_no_crash)
{
    MonomerToolTestFixture fix;

    // place a monomer
    auto monomer_pos = QPointF(100, 100);
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    fix.mouseClick(monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // hover over the monomer so the tool creates unbound AP items
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    fix.mouseMove(monomer_pos);

    // hover over an attachment point so the tool also draws a hint fragment
    auto c_ap_pos = fix.getAttachmentPointPos(0, "C");
    fix.mouseMove(c_ap_pos);

    // undo the placement (simulates Ctrl+Z). This destroys the monomer
    // graphics item and, transitively, its child AP items.
    auto* undo_stack = fix.m_mol_model->findChild<QUndoStack*>();
    BOOST_REQUIRE(undo_stack != nullptr);
    undo_stack->undo();
    process_qt_events();
    fix.confirmIsEmpty();

    // moving the mouse a little crashes pre-fix because the tool's cached
    // m_unbound_ap_items / m_hovered_monomeric_item are dangling
    fix.mouseMove(c_ap_pos + QPointF(5, 5));
}

/**
 * Start a drag from one monomer toward another so the tool caches drag-end AP
 * items parented to the second monomer, then undo (which destroys the second
 * monomer) before releasing. The tool's onStructureUpdated() must drop its
 * stale drag-end pointers without trying to delete them - the items were
 * already destroyed as Qt children of the now-gone monomer.
 *
 * NOTE: catches the bug reliably only under ASan; a plain debug build may
 * pass even with the bug present if the freed memory hasn't been reused.
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

    // start a drag from A's C attachment point
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    fix.mouseMove(pos_a);
    auto a_c_ap = fix.getAttachmentPointPos(0, "C");
    fix.mouseMove(a_c_ap);
    fix.mousePress(a_c_ap);

    // drag over B - this populates m_drag_end_unbound_ap_items, parented to B
    fix.mouseMove(pos_b, Qt::LeftButton);
    auto b_n_ap = fix.getAttachmentPointPos(1, "N");
    fix.mouseMove(b_n_ap, Qt::LeftButton);

    // undo the second click - destroys monomer B and (transitively) the
    // drag-end AP items parented to it. The structure-update callback runs
    // synchronously here; pre-fix-extension this hits a double-free inside
    // clearDragEndAttachmentPointsLabels().
    auto* undo_stack = fix.m_mol_model->findChild<QUndoStack*>();
    BOOST_REQUIRE(undo_stack != nullptr);
    undo_stack->undo();
    process_qt_events();
    BOOST_TEST(fix.m_mol_model->getMol()->getNumAtoms() == 1);

    // releasing the drag must not crash either
    fix.mouseRelease(b_n_ap);
}

/**
 * Make sure that dragging from a nucleic acid base's "pair" attachment point to
 * a peptide monomer doesn't cause a crash.
 */
BOOST_AUTO_TEST_CASE(test_drag_from_na_base_to_peptide)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}|RNA1{A}$$$$V2.0");
    fix.setMonomericConnectionTool(
        MonomericConnectionTool::COVALENT_OR_DISULFIDE);
    auto start_monomer_pos = fix.getMonomerPos(1);
    auto end_monomer_pos = fix.getMonomerPos(0);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(1, "H-bond");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    fix.mouseRelease(end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}|RNA1{A}$RNA1,PEPTIDE1,1:pair-1:pair$$$V2.0");
}

/**
 * Make sure that dragging to a nucleic acid base from a peptide monomer doesn't
 * cause a crash (i.e. dragging in the opposite direction relative to the last
 * test).
 */
BOOST_AUTO_TEST_CASE(test_drag_to_na_base_from_peptide)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}|RNA1{A}$$$$V2.0");
    fix.setNucleicAcidTool(NucleicAcidTool::A);
    auto start_monomer_pos = fix.getMonomerPos(0);
    auto end_monomer_pos = fix.getMonomerPos(1);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(0, "X");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    fix.mouseRelease(end_monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}|RNA1{A}$PEPTIDE1,RNA1,1:pair-1:pair$$$V2.0");
}

} // namespace sketcher
} // namespace schrodinger
