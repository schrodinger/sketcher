#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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
    TestSketcherWidget sk;
    auto mol_model = sk.m_mol_model;
    TestBackgroundContextMenu menu(sk.m_sketcher_model);

    menu.updateActions();
    BOOST_TEST(!menu.m_undo_act->isEnabled());
    BOOST_TEST(!menu.m_redo_act->isEnabled());
    check_nonempty_actions(menu, false);

    import_mol_text(mol_model, "CCCC");
    // Use updateActions() as a placeholder for showEvent()
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
