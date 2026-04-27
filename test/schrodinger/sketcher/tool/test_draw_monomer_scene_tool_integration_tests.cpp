#define BOOST_TEST_MODULE test_draw_monomer_scene_tool_integration_tests

#include <memory>
#include <string>
#include <vector>

#include <QApplication>
#include <QEvent>
#include <QGraphicsSceneMouseEvent>
#include <QPointF>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Process all Qt events, includeing DeferredDelete events.
 */
void process_qt_events()
{
    // call processEvents multiple times in case an any current events put new
    // events on the queue (e.g. starting a timer with a timeout of 0)
    for (int i = 0; i < 3; ++i) {
        QApplication::processEvents();
        // processEvents will never process DeferredDelete events, but
        // sendPostedEvents will if we explicitly pass their event type.
        // (Despite what Qt's documentation claims, passing 0 as the event type
        // processes everything *other than* DeferredDeletes.)
        QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    }
}

/**
 * Set both the scene and screen pos for an event
 */
void set_event_pos(QGraphicsSceneMouseEvent& event, const QPointF& pos)
{
    event.setScenePos(pos);
    event.setScreenPos(pos.toPoint());
}

/**
 * Fixture that provides a clean Scene with MolModel and SketcherModel for
 * testing monomer drawing tools.
 */
struct MonomerToolTestFixture {
    std::shared_ptr<TestScene> m_scene;
    MolModel* m_mol_model;
    SketcherModel* m_sketcher_model;

    MonomerToolTestFixture()
    {
        m_scene = TestScene::getScene();
        m_mol_model = m_scene->m_mol_model;
        m_sketcher_model = m_scene->m_sketcher_model;
        m_sketcher_model->setValue(
            ModelKey::INTERFACE_TYPE,
            static_cast<int>(InterfaceType::ATOMISTIC_OR_MONOMERIC));
        process_qt_events();
    }

    void setAminoAcidTool(AminoAcidTool tool)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::AMINO_ACID)},
             {ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(tool)},
             {ModelKey::AMINO_ACID_SYMBOL, QString("")}});
        process_qt_events();
    }

    void setNucleicAcidTool(NucleicAcidTool tool)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
             {ModelKey::NUCLEIC_ACID_TOOL, QVariant::fromValue(tool)}});
        process_qt_events();
    }

    void importMolText(const std::string& text)
    {
        import_mol_text(m_mol_model, text);
        process_qt_events();
    }

    QPointF getMonomerPos(unsigned int monomer_idx)
    {
        auto mol = m_mol_model->getMol();
        BOOST_REQUIRE(mol);
        BOOST_REQUIRE(monomer_idx < mol->getNumAtoms());
        return to_scene_xy(mol->getConformer().getAtomPos(monomer_idx));
    }

    QPointF getAttachmentPointPos(unsigned int monomer_idx,
                                  const std::string& ap_display_name)
    {
        auto monomer_pos = getMonomerPos(monomer_idx);
        auto* monomer_item = m_scene->getTopInteractiveItemAt(
            monomer_pos, InteractiveItemFlag::MONOMER);
        BOOST_REQUIRE(monomer_item != nullptr);

        // Search through child items to find the attachment point
        for (auto* child : monomer_item->childItems()) {
            auto* ap_item =
                qgraphicsitem_cast<UnboundMonomericAttachmentPointItem*>(child);
            if (ap_item) {
                auto ap = ap_item->getAttachmentPoint();
                if (ap.display_name == ap_display_name) {
                    auto offset_vec = direction_to_qt_vector(ap.direction);
                    auto offset_dist = UNBOUND_AP_LINE_LENGTH - 1 +
                                       STANDARD_AA_BORDER_WIDTH / 2;
                    return monomer_pos + offset_dist * offset_vec;
                }
            }
        }
        throw std::runtime_error("Attachment point " + ap_display_name +
                                 " not found");
    }

    void verifyHELM(const std::string& expected)
    {
        auto actual = get_mol_text(m_mol_model, rdkit_extensions::Format::HELM);
        BOOST_TEST(actual == expected);
    }

    void mouseMove(const QPointF& pos,
                   const Qt::MouseButtons btns = Qt::NoButton)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMouseMove);
        set_event_pos(event, pos);
        event.setButton(Qt::NoButton);
        event.setButtons(btns);
        m_scene->mouseMoveEvent(&event);
        process_qt_events();
    }

    void mousePress(const QPointF& pos)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMousePress);
        set_event_pos(event, pos);
        event.setButton(Qt::LeftButton);
        event.setButtons(Qt::LeftButton);
        m_scene->mousePressEvent(&event);
        process_qt_events();
    }

    void mouseRelease(const QPointF& pos)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMouseRelease);
        set_event_pos(event, pos);
        event.setButton(Qt::LeftButton);
        event.setButtons(Qt::NoButton);
        m_scene->mouseReleaseEvent(&event);
        process_qt_events();
    }

    void mouseClick(const QPointF& pos)
    {
        mouseMove(pos);
        mousePress(pos);
        mouseRelease(pos);
    }

    void mouseDrag(const QPointF& start, const QPointF& end)
    {
        mouseMove(start);
        mousePress(start);
        mouseMove(end, Qt::LeftButton);
        mouseRelease(end);
    }

    void confirmIsEmpty()
    {
        auto mol = m_mol_model->getMol();
        BOOST_TEST(mol->getNumAtoms() == 0);
    }
};

/**
 * Confirm that clicking in empty space adds the appropriate monomer
 */
BOOST_AUTO_TEST_CASE(test_click_empty_space_adds_monomer)
{
    MonomerToolTestFixture fix;
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    fix.mouseClick({0, 0});
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that clicking on an existing monomer with the equivalent monomer tool
 * adds a new monomer via the default attachment point.
 */
BOOST_AUTO_TEST_CASE(test_click_existing_monomer_same_residue_adds_residue)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto pos = fix.getMonomerPos(0);
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    fix.mouseClick(pos);
    fix.verifyHELM("PEPTIDE1{A.A}$$$$V2.0");
}

/**
 * Confirm that clicking on an existing monomer with a different monomer tool of
 * the same monomer type mutates the monomer.
 */
BOOST_AUTO_TEST_CASE(test_click_existing_monomer_different_residue_mutates)
{
    MonomerToolTestFixture fix;

    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto pos = fix.getMonomerPos(0);
    fix.setAminoAcidTool(AminoAcidTool::CYS);
    fix.mouseClick(pos);
    fix.verifyHELM("PEPTIDE1{C}$$$$V2.0");
}

/**
 * Confirm that clicking on an existing monomer with a monomer tool of
 * a different monomer type has no effect.
 */
BOOST_AUTO_TEST_CASE(test_click_existing_monomer_different_monomer_type)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto pos = fix.getMonomerPos(0);
    fix.setNucleicAcidTool(NucleicAcidTool::P);
    fix.mouseClick(pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");
}

/**
 * Confirm that clicking on an unbounnd attachment point of an existing monomer
 * adds a new monomer via the clicked attachment point.
 */
BOOST_AUTO_TEST_CASE(test_click_attachment_point)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    auto monomer_pos = fix.getMonomerPos(0);
    fix.setAminoAcidTool(AminoAcidTool::CYS);
    // hover over the monomer to trigger AP label creation
    fix.mouseMove(monomer_pos);

    // click on the N terminus attachment point
    fix.setAminoAcidTool(AminoAcidTool::CYS);
    fix.mouseMove(monomer_pos);
    auto n_ap_pos = fix.getAttachmentPointPos(0, "N");
    fix.mouseClick(n_ap_pos);
    fix.verifyHELM("PEPTIDE1{C.A}$$$$V2.0");

    // click on the C terminus attachment point
    fix.setAminoAcidTool(AminoAcidTool::PHE);
    fix.mouseMove(monomer_pos);
    auto c_ap_pos = fix.getAttachmentPointPos(0, "C");
    fix.mouseClick(c_ap_pos);
    fix.verifyHELM("PEPTIDE1{C.A.F}$$$$V2.0");

    // click on the side chain attachment point
    fix.setAminoAcidTool(AminoAcidTool::TRP);
    fix.mouseMove(monomer_pos);
    auto x_ap_pos = fix.getAttachmentPointPos(0, "X");
    fix.mouseClick(x_ap_pos);
    fix.verifyHELM(
        "PEPTIDE1{C.A.F}|PEPTIDE2{W}$PEPTIDE1,PEPTIDE2,2:R3-1:R3$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer adds a new monomer using
 * the default attachment point
 */
BOOST_AUTO_TEST_CASE(test_drag_monomer_to_empty_adds_connected_default_ap)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    fix.setAminoAcidTool(AminoAcidTool::CYS);

    // Drag from the monomer to the right
    auto start_pos = fix.getMonomerPos(0);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A.C}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer adds a new monomer using
 * the default attachment point, even when the monomer tool is equivalent to the
 * existing monomer (e.g. ALA tool on an A monomer).
 */
BOOST_AUTO_TEST_CASE(
    test_drag_monomer_to_empty_adds_connected_default_ap_same_monomer)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{A}$$$$V2.0");
    fix.setAminoAcidTool(AminoAcidTool::ALA);

    // Drag from the monomer to the right
    auto start_pos = fix.getMonomerPos(0);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{A.A}$$$$V2.0");
}

/**
 * Confirm that click-and-drag from the attachment point of an existing monomer
 * to empty space adds a new monomer via the specified attachment point of the
 * existing monomer
 */
BOOST_AUTO_TEST_CASE(test_drag_ap_to_empty_adds_connected_via_dragged_ap)
{
    MonomerToolTestFixture fix;
    fix.setAminoAcidTool(AminoAcidTool::CYS);

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

    fix.verifyHELM("PEPTIDE1{C.A}$$$$V2.0");

    // hover over the first monomer so that its attachment point graphics items
    // are created again
    fix.setAminoAcidTool(AminoAcidTool::PHE);
    fix.mouseMove(ala_pos);

    // Drag from C attachment point to empty space
    start_pos = fix.getAttachmentPointPos(0, "C");
    end_pos = start_pos + QPointF(100, 100);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM("PEPTIDE1{C.A.F}$$$$V2.0");

    // hover over the first monomer so that its attachment point graphics items
    // are created again
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    fix.mouseMove(ala_pos);

    // Drag from C attachment point to empty space
    start_pos = fix.getAttachmentPointPos(0, "X");
    end_pos = start_pos + QPointF(-50, 100);
    fix.mouseDrag(start_pos, end_pos);

    fix.verifyHELM(
        "PEPTIDE1{C.A.F}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,2:R3-1:R3$$$V2.0");
}

/**
 * Confirm that click-and-drag from an existing monomer to an existing monomer
 * connects them via the default attachment points
 */
BOOST_AUTO_TEST_CASE(test_drag_monomer_to_monomer_connects_default_aps)
{
    MonomerToolTestFixture fix;
    fix.setAminoAcidTool(AminoAcidTool::ALA);
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
    fix.setAminoAcidTool(AminoAcidTool::ALA);
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
 * Confirm that dragging from empty space to empty space with a peptide monomer
 * creates a dimer with a backbone connection
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_empty_adds_two_connected_default_aps)
{
    MonomerToolTestFixture fix;
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    auto start_pos = QPointF(100, 100);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);
    fix.verifyHELM("PEPTIDE1{A.A}$$$$V2.0");
}

/**
 * Confirm that dragging from empty space to empty space with a nucleic acid
 * base monomer creates two paired bases
 */
BOOST_AUTO_TEST_CASE(test_nucleic_acid_base_drag_empty_to_empty_uses_pair_ap)
{
    MonomerToolTestFixture fix;
    fix.setNucleicAcidTool(NucleicAcidTool::A);
    auto start_pos = QPointF(100, 100);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);
    fix.verifyHELM("RNA1{A}|RNA2{A}$RNA1,RNA2,1:pair-1:pair$$$V2.0");
}

/**
 * Confirm that dragging from empty space to empty space with a nucleic acid
 * phosphate monomer is ignored and doesn't create any monomers
 */
BOOST_AUTO_TEST_CASE(test_nucleic_acid_sugar_drag_empty_to_empty_ignored)
{
    MonomerToolTestFixture fix;
    fix.setNucleicAcidTool(NucleicAcidTool::R);
    auto start_pos = QPointF(100, 100);
    auto end_pos = start_pos + QPointF(100, 0);
    fix.mouseDrag(start_pos, end_pos);
    fix.confirmIsEmpty();
}

/**
 * Confirm that dragging from empty space to an existing monomer creates a new
 * monomer and connects it to the default attachment point of the existing
 * monomer
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_monomer_adds_connected_default_aps)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{C}$$$$V2.0");
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    auto start_pos = QPointF(100, 100);
    auto end_pos = fix.getMonomerPos(0);
    fix.mouseDrag(start_pos, end_pos);
    fix.verifyHELM("PEPTIDE1{A.C}$$$$V2.0");
}

/**
 * Confirm that dragging from empty space to an existing monomer creates a new
 * monomer and connects it to the specified attachment point of the existing
 * monomer
 */
BOOST_AUTO_TEST_CASE(test_drag_empty_to_ap_adds_connected_correct_aps)
{
    MonomerToolTestFixture fix;
    fix.importMolText("PEPTIDE1{C}$$$$V2.0");
    fix.setAminoAcidTool(AminoAcidTool::ALA);
    auto start_pos = QPointF(100, 100);
    auto monomer_pos = fix.getMonomerPos(0);
    fix.mouseMove(start_pos);
    fix.mousePress(start_pos);
    // first, drag to the existing monomer to make its attachment points appear
    fix.mouseMove(monomer_pos);
    auto end_pos = fix.getAttachmentPointPos(0, "X");
    fix.mouseMove(end_pos);
    fix.mouseRelease(end_pos);
    fix.verifyHELM(
        "PEPTIDE1{C}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,1:R3-1:R2$$$V2.0");
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
    fix.setAminoAcidTool(AminoAcidTool::CYS);
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
    // the drag should have added a new monomer to the N terminus of the start
    // monomer since we released over an invalid attachment point
    fix.verifyHELM("PEPTIDE1{C.A.A}$$$$V2.0");
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
    fix.setAminoAcidTool(AminoAcidTool::PHE);
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
    // the drag should have added a new monomer with a side chain interaction
    // since we released over an invalid attachment point
    fix.verifyHELM("PEPTIDE1{C}|PEPTIDE2{A}|PEPTIDE3{F}$PEPTIDE1,PEPTIDE2,1:R2-"
                   "1:R2|PEPTIDE1,PEPTIDE3,1:R3-1:R3$$$V2.0");
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
    fix.setAminoAcidTool(AminoAcidTool::CYS);
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
    // the drag should have added a new monomer with a side chain interaction
    // since we released over an invalid attachment point
    fix.verifyHELM("PEPTIDE1{C.C.C}$PEPTIDE1,PEPTIDE1,2:R3-3:R3$$$V2.0");
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
    fix.setAminoAcidTool(AminoAcidTool::ALA);

    // place a monomer
    auto monomer_pos = QPointF(100, 100);
    fix.mouseClick(monomer_pos);
    fix.verifyHELM("PEPTIDE1{A}$$$$V2.0");

    // hover over the monomer so the tool creates unbound AP items
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
    // m_unbound_ap_items / m_hovered_item are dangling
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
    fix.setAminoAcidTool(AminoAcidTool::ALA);

    // place two monomers far enough apart that they don't overlap
    fix.mouseClick(QPointF(100, 100));
    auto pos_a = fix.getMonomerPos(0);
    fix.mouseClick(pos_a + QPointF(200, 0));
    BOOST_REQUIRE(fix.m_mol_model->getMol()->getNumAtoms() == 2);
    auto pos_b = fix.getMonomerPos(1);

    // start a drag from A's C attachment point
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
    fix.setNucleicAcidTool(NucleicAcidTool::A);
    auto start_monomer_pos = fix.getMonomerPos(1);
    auto end_monomer_pos = fix.getMonomerPos(0);
    fix.mouseMove(start_monomer_pos);
    auto start_ap_pos = fix.getAttachmentPointPos(1, "pair");
    fix.mouseMove(start_ap_pos);
    fix.mousePress(start_ap_pos);
    fix.mouseMove(end_monomer_pos);
    fix.mouseRelease(end_monomer_pos);
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
