#define BOOST_TEST_MODULE Test_Sketcher
#include <QAction>
#include <QButtonGroup>
#include <QKeyEvent>
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_top_bar.h"
#include "schrodinger/sketcher/widget/sketcher_top_bar.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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
    TestSketcherTopBar(SketcherModel* model)
    {
        setModel(model);
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
        bool exp_value = model->getValueBool(signal_packet.key);
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
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestSketcherTopBar top_bar(model);
    import_mol_text(test_scene->m_mol_model, "CC");

    // We must call this method manually here because there is no signal from
    // the sketcher scene during unit tests
    top_bar.updateWidgetsEnabled();

    // setModel() should synchronize the view and model immediately
    BOOST_TEST(view_synchronized_to_model(top_bar));

    // Try interacting with subwidgets of the view. Synchronization should occur
    // automatically.
    for (auto& signal_packet : top_bar.getSignalPackets()) {
        bool value = model->getValueBool(signal_packet.key);
        signal_packet.action->trigger();
        BOOST_TEST(signal_packet.action->isChecked() != value);
        BOOST_TEST(view_synchronized_to_model(top_bar));
    }

    // Try changing the state of the model. The view should update the state of
    // its subwidgets automatically.
    for (auto& setter_packet : top_bar.getSetterPackets()) {
        auto new_value = !model->getValueBool(setter_packet.key);
        model->setValue(setter_packet.key, new_value);
        BOOST_TEST(model->getValueBool(setter_packet.key) == new_value);
        BOOST_TEST(view_synchronized_to_model(top_bar));
    }
}

/**
 * Verify that scene-dependent QActions/QWidgets are enabled and disabled
 * appropriately in response to changes in the model.
 */
BOOST_AUTO_TEST_CASE(widgets_enabled)
{
    TestSketcherWidget sk;
    auto model = sk.m_sketcher_model;
    auto mol_model = sk.m_mol_model;
    TestSketcherTopBar top_bar(model);
    auto menu = top_bar.getMoreActionsMenu();
    import_mol_text(mol_model, "CC");

    mol_model->selectAll();
    auto atoms = mol_model->getSelectedAtoms();

    for (auto sel_empty : {true, false}) {
        if (sel_empty) {
            mol_model->clearSelection();
        } else {
            mol_model->select({atoms}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
        }
        BOOST_TEST(menu->m_clear_selection_act->isEnabled() == !sel_empty);
        for (auto has_rxn : {true, false}) {
            if (has_rxn) {
                mol_model->addNonMolecularObject(
                    NonMolecularType::RXN_ARROW,
                    RDGeom::Point3D(0.0, 0.0, 0.0));
            }
            BOOST_TEST(model->hasReaction() == has_rxn);
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
                sk.m_undo_stack->undo();
            }
            BOOST_TEST(!model->hasReaction());
        }
    }
}

BOOST_AUTO_TEST_CASE(testhandleShortcutAction)
{
    TestSketcherWidget sk;
    TestSketcherTopBar top_bar(sk.m_sketcher_model);
    auto mol_model = sk.m_mol_model;
    // Populate the undo and redo stacks so that those buttons will be enabled
    import_mol_text(mol_model, "CC");
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

    mol_model->selectAll();
    BOOST_CHECK_EQUAL(menu->m_select_all_act->isEnabled(), false);
    BOOST_CHECK_EQUAL(menu->m_clear_selection_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(clear_all_seq), true);
    BOOST_CHECK_EQUAL(clear_spy.count(), 1);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(select_all_seq), true);
    BOOST_CHECK_EQUAL(select_all_spy.count(), 1);

    QKeySequence undo_seq{Qt::CTRL | Qt::Key_Z};
    QSignalSpy undo_spy{&top_bar, &SketcherTopBar::undoRequested};
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(undo_seq), true);
    BOOST_CHECK_EQUAL(undo_spy.count(), 1);

    QKeySequence redo_seq{Qt::CTRL | Qt::SHIFT | Qt::Key_Z};
    QSignalSpy redo_spy{&top_bar, &SketcherTopBar::redoRequested};
    BOOST_CHECK_EQUAL(menu->m_redo_act->isEnabled(), false);
    sk.m_undo_stack->undo();
    BOOST_CHECK_EQUAL(menu->m_redo_act->isEnabled(), true);
    BOOST_CHECK_EQUAL(top_bar.handleShortcutAction(redo_seq), true);
    BOOST_CHECK_EQUAL(redo_spy.count(), 1);
}

} // namespace sketcher
} // namespace schrodinger
