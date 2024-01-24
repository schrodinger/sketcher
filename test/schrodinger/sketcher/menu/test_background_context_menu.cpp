#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "../test_sketcherScene.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

class TestBackgroundContextMenu : public BackgroundContextMenu
{
  public:
    TestBackgroundContextMenu(SketcherModel* model) :
        BackgroundContextMenu(model)
    {
    }
    using BackgroundContextMenu::m_cut_copy_manager;
    using BackgroundContextMenu::m_export_to_file_act;
    using BackgroundContextMenu::m_redo_act;
    using BackgroundContextMenu::m_save_image_act;
    using BackgroundContextMenu::m_select_all_act;
    using BackgroundContextMenu::m_undo_act;
    using BackgroundContextMenu::updateActions;
};

void check_nonempty_actions(TestBackgroundContextMenu& menu, bool exp_enabled)
{
    BOOST_TEST(menu.m_save_image_act->isEnabled() == exp_enabled);
    BOOST_TEST(menu.m_export_to_file_act->isEnabled() == exp_enabled);
    BOOST_TEST(menu.m_select_all_act->isEnabled() == exp_enabled);
    BOOST_TEST(menu.m_cut_copy_manager->m_copy_action->isEnabled() ==
               exp_enabled);
    BOOST_TEST(menu.m_cut_copy_manager->m_copy_as_menu->isEnabled() ==
               exp_enabled);
}

/**
 * Verify that the menu is properly updated to match the state of the model.
 */
BOOST_AUTO_TEST_CASE(test_updateActions_deprecated)
{
    testSketcherScene scene;
    auto model = scene.getModel();
    TestBackgroundContextMenu menu(model);
    scene.connectContextMenu(menu);

    menu.updateActions();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    scene.importText("CCCC");
    // Use updateActions() as a placeholder for showEvent()
    menu.updateActions();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    scene.undo();
    menu.updateActions();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    scene.redo();
    menu.updateActions();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    BOOST_TEST(menu.m_select_all_act->isEnabled() == true);
}

/**
 * Verify that the menu is properly updated to match the state of the model.
 */
BOOST_AUTO_TEST_CASE(test_updateActions)
{
    TestSketcherWidget sk;
    TestBackgroundContextMenu menu(sk.m_sketcher_model);

    menu.updateActions();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    sk.addFromString("CCCC");
    menu.updateActions();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    sk.m_undo_stack->undo();
    menu.updateActions();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    sk.m_undo_stack->redo();
    menu.updateActions();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    BOOST_TEST(menu.m_select_all_act->isEnabled() == true);
}

} // namespace sketcher
} // namespace schrodinger
