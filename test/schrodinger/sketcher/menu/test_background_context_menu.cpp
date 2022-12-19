#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>
#include "../test_common.h"
#include "../test_sketcherScene.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/cut_copy_action_manager.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"

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
BOOST_AUTO_TEST_CASE(test_updateActions)
{
    testSketcherScene scene;
    auto model = scene.getModel();
    TestBackgroundContextMenu menu(model);
    scene.connectContextMenu(menu);

    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    scene.importText("CCCC");
    emit model->sceneContentsChanged();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    scene.undo();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    scene.redo();
    BOOST_TEST(menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, true);

    BOOST_TEST(menu.m_select_all_act->isEnabled() == true);
}

} // namespace sketcher
} // namespace schrodinger
