#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Verify that cut and copy actions are enabled as expected
 */
BOOST_AUTO_TEST_CASE(test_updateActions)
{
    SketcherModel model;
    CutCopyActionManager mgr(nullptr);
    sketcherScene scene;
    mgr.setModel(&model);
    scene.setModel(&model);
    // empty
    BOOST_TEST(!mgr.m_cut_action->isEnabled());
    BOOST_TEST(!mgr.m_copy_action->isEnabled());
    BOOST_TEST(!mgr.m_copy_as_menu->isEnabled());

    scene.importText("CC");
    BOOST_TEST(!mgr.m_cut_action->isEnabled());
    BOOST_TEST(mgr.m_copy_action->isEnabled());
    BOOST_TEST(mgr.m_copy_as_menu->isEnabled());
    BOOST_TEST(mgr.m_copy_action->text().toStdString() == "Copy All");
    BOOST_TEST(mgr.m_copy_as_menu->title().toStdString() == "Copy All As");

    scene.selectAll();
    BOOST_TEST(mgr.m_cut_action->isEnabled());
    BOOST_TEST(mgr.m_copy_action->isEnabled());
    BOOST_TEST(mgr.m_copy_as_menu->isEnabled());
    BOOST_TEST(mgr.m_copy_action->text().toStdString() == "Copy");
    BOOST_TEST(mgr.m_copy_as_menu->title().toStdString() == "Copy As");

    // All formats are present as actions, we just show/hide based on whether
    // there is a mol or a reaction present
    BOOST_TEST(mgr.m_copy_as_menu->actions().size() ==
               get_standard_export_formats().size() +
                   get_reaction_export_formats().size());

    // confirm copy as menu toggles based on reactions
    auto reaction_actions_visible = [&mgr](bool expect_reaction) {
        for (auto act : mgr.m_copy_as_menu->actions()) {
            bool expected = act->data().toBool() == expect_reaction;
            if (act->isVisible() != expected) {
                return false;
            }
        }
        return true;
    };

    auto select_one_atom = [&scene]() {
        auto atoms = scene.quickGetAtoms();
        std::unordered_set<sketcherGraphicalObject*> sel_objs = {atoms.front()};
        scene.setSelection(sel_objs);
    };

    BOOST_TEST(reaction_actions_visible(false));
    scene.clearStructure();
    BOOST_TEST(reaction_actions_visible(false));
    scene.importText("CC>>CC");
    BOOST_TEST(reaction_actions_visible(true));
    select_one_atom();
    BOOST_TEST(reaction_actions_visible(false));
    scene.selectAll();
    BOOST_TEST(reaction_actions_visible(true));
    scene.clearStructure();
    BOOST_TEST(reaction_actions_visible(false));
    scene.importText("CC");
    BOOST_TEST(reaction_actions_visible(false));

    // background context menu always "Copy All"
    mgr.setAlwaysCopyAll(true);

    BOOST_TEST(reaction_actions_visible(false));
    scene.clearStructure();
    BOOST_TEST(reaction_actions_visible(false));
    scene.importText("CC>>CC");
    BOOST_TEST(reaction_actions_visible(true));
    select_one_atom();
    BOOST_TEST(reaction_actions_visible(true)); // ignores selection
    scene.selectAll();
    BOOST_TEST(reaction_actions_visible(true)); // ignores selection
    scene.clearStructure();
    BOOST_TEST(reaction_actions_visible(false));
    scene.importText("CC");
    BOOST_TEST(reaction_actions_visible(false));
}

} // namespace sketcher
} // namespace schrodinger
