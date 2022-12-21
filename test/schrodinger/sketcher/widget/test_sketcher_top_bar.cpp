#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>
#include "../test_common.h"
#include "schrodinger/sketcher/widget/sketcher_top_bar.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_top_bar.h"
#include "schrodinger/sketcher/cut_copy_action_manager.h"
#include "schrodinger/sketcher/reaction_arrow.h"
#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"
#include "schrodinger/sketcher/Scene.h"
#include "../test_sketcherScene.h"
#include <QSignalSpy>
#include <QTest>
#include <QKeyEvent>
#include <QButtonGroup>
#include <QAction>

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{
std::vector<SetterPacket> m_setter_packets;

/**
 * Data necessary for updating the state of the model from checkable
 * actions on this view.
 *
 * These are used when synchronizing the model to this view.
 */
std::vector<SignalPacket> m_signal_packets;

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestSketcherTopBar : public SketcherTopBar
{
  public:
    TestSketcherTopBar()
    {
        setModel(new SketcherModel(this));
    }

    std::vector<SignalPacket> getSignalPackets() const
    {
        return m_signal_packets;
    }

    ImportMenu* getImportMenu() const
    {
        return m_import_menu;
    }

    MoreActionsMenu* getMoreActionsMenu() const
    {
        return m_more_actions_menu;
    }
    using SketcherTopBar::updateWidgetsEnabled;
};

/**
 * Compare the state of the top bar's subwidgets and its corresponding model.
 *
 * @param top_bar A top bar with an assigned model
 */
bool view_synchronized_to_model(TestSketcherTopBar& top_bar)
{
    auto model = top_bar.getModel();
    for (auto& signal_packet : top_bar.getSignalPackets()) {
        bool exp_value = model->getValue(signal_packet.key).toBool();
        if (signal_packet.action->isChecked() != exp_value) {
            return false;
        }
    }
    return true;
}

/**
 * Verify proper synchronization with model.
 */
BOOST_AUTO_TEST_CASE(synchronize)
{

    TestSketcherTopBar top_bar;
    auto model = top_bar.getModel();
    sketcherScene scene;
    scene.setModel(model);
    scene.importText("CC");

    // We must call this method manually here because there is no signal from
    // the sketcher scene during unit tests
    top_bar.updateWidgetsEnabled();

    // setModel() should synchronize the view and model immediately
    BOOST_TEST(view_synchronized_to_model(top_bar));

    // Try interacting with subwidgets of the view. Synchronization should occur
    // automatically.
    for (auto& signal_packet : top_bar.getSignalPackets()) {
        bool value = model->getValue(signal_packet.key).toBool();
        signal_packet.action->trigger();
        BOOST_TEST(signal_packet.action->isChecked() != value);
        BOOST_TEST(view_synchronized_to_model(top_bar));
    }

    // Try changing the state of the model. The view should update the state of
    // its subwidgets automatically.
    for (auto& setter_packet : top_bar.getSetterPackets()) {
        int new_value =
            static_cast<int>(!model->getValue(setter_packet.key).toInt());
        model->setValue(setter_packet.key, QVariant(new_value));
        BOOST_TEST(model->getValue(setter_packet.key).toInt() == new_value);
        BOOST_TEST(view_synchronized_to_model(top_bar));
    }
}

/**
 * Verify that scene-dependent QActions/QWidgets are enabled and disabled
 * appropriately in response to changes in the model.
 */
BOOST_AUTO_TEST_CASE(widgets_enabled)
{

    TestSketcherTopBar top_bar;
    auto model = top_bar.getModel();
    auto menu = top_bar.getMoreActionsMenu();
    sketcherScene scene;
    scene.setModel(model);
    scene.importText("CC");

    auto atoms = scene.quickGetAtoms();
    std::unordered_set<QGraphicsItem*> atom_sel(atoms.begin(), atoms.end());

    QGraphicsTextItem text_item;
    std::unordered_set<QGraphicsItem*> selection = {&text_item};

    for (auto sel_empty : {true, false}) {
        if (sel_empty) {
            scene.clearSelection();
        } else {
            scene.setSelection(atom_sel);
        }
        BOOST_TEST(menu->m_clear_selection_act->isEnabled() == !sel_empty);
        for (auto has_rxn : {true, false}) {
            if (has_rxn) {
                scene.addReactionArrowAt(QPointF(0, 0));
            }
            BOOST_TEST(scene.numberOfReactionSteps() == int(has_rxn));
            auto cc = menu->m_cut_copy_manager;
            BOOST_TEST(cc->m_cut_action->isEnabled() == !sel_empty);
            BOOST_TEST(cc->m_copy_action->isEnabled() == true);
            BOOST_TEST(cc->m_copy_as_menu->menuAction()->isEnabled() == true);
            auto copy_text = cc->m_copy_action->text().toStdString();
            auto copy_as_text =
                cc->m_copy_as_menu->menuAction()->text().toStdString();
            if (sel_empty) {
                BOOST_TEST(copy_text == "Copy All");
                BOOST_TEST(copy_as_text == "Copy All As");
            } else {
                BOOST_TEST(copy_text == "Copy");
                BOOST_TEST(copy_as_text == "Copy As");
            }

            if (has_rxn) {
                // Remove the reaction arrow for further testing
                scene.undo();
            }
            BOOST_TEST(scene.numberOfReactionSteps() == 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(testhandleShortcutAction)
{
    TestSketcherTopBar top_bar;
    sketcherScene scene;
    scene.setModel(top_bar.getModel());
    // Populate the undo and redo stacks so that those buttons will be enabled
    scene.importText("CC");
    top_bar.updateWidgetsEnabled();

    QKeySequence select_all_seq{Qt::CTRL | Qt::Key_A};
    QKeySequence clear_all_seq{Qt::CTRL | Qt::Key_D};
    QSignalSpy select_all_spy{&top_bar, &SketcherTopBar::selectAllRequested};
    QSignalSpy clear_spy{&top_bar, &SketcherTopBar::clearSelectionRequested};
    auto menu = top_bar.getMoreActionsMenu();

    BOOST_CHECK_EQUAL(menu->m_select_all_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(menu->m_clear_selection_act->isEnabled(), false);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(clear_all_seq), true);
    BOOST_CHECK_EQUAL(clear_spy.count(), 0);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(select_all_seq), true);
    BOOST_CHECK_EQUAL(select_all_spy.count(), 1);

    scene.selectAll();
    BOOST_CHECK_EQUAL(menu->m_select_all_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(menu->m_clear_selection_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(clear_all_seq), true);
    BOOST_CHECK_EQUAL(clear_spy.count(), 1);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(select_all_seq), true);
    BOOST_CHECK_EQUAL(select_all_spy.count(), 2);

    QKeySequence undo_seq{Qt::CTRL | Qt::Key_Z};
    QSignalSpy undo_spy{&top_bar, &SketcherTopBar::undoRequested};
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(undo_seq), true);
    BOOST_CHECK_EQUAL(undo_spy.count(), 1);

    QKeySequence redo_seq{Qt::CTRL | Qt::SHIFT | Qt::Key_Z};
    QSignalSpy redo_spy{&top_bar, &SketcherTopBar::redoRequested};
    BOOST_CHECK_EQUAL(menu->m_redo_act->isEnabled(), false);
    scene.undo();
    BOOST_CHECK_EQUAL(menu->m_redo_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(redo_seq), true);
    BOOST_CHECK_EQUAL(redo_spy.count(), 1);
}

} // namespace sketcher
} // namespace schrodinger
