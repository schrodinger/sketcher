#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/mol_model.h"
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
    TestSketcherWidget sk;
    auto mol_model = sk.m_mol_model;
    CutCopyActionManager mgr(nullptr);
    mgr.setModel(sk.m_sketcher_model);
    // empty
    BOOST_TEST(!mgr.m_cut_action->isEnabled());
    BOOST_TEST(!mgr.m_copy_action->isEnabled());
    BOOST_TEST(!mgr.m_copy_as_menu->isEnabled());

    sk.addFromString("CC");
    BOOST_TEST(!mgr.m_cut_action->isEnabled());
    BOOST_TEST(mgr.m_copy_action->isEnabled());
    BOOST_TEST(mgr.m_copy_as_menu->isEnabled());
    BOOST_TEST(mgr.m_copy_action->text().toStdString() == "Copy All");
    BOOST_TEST(mgr.m_copy_as_menu->title().toStdString() == "Copy All As");

    mol_model->selectAll();
    BOOST_TEST(mgr.m_cut_action->isEnabled());
    BOOST_TEST(mgr.m_copy_action->isEnabled());
    BOOST_TEST(mgr.m_copy_as_menu->isEnabled());
    BOOST_TEST(mgr.m_copy_action->text().toStdString() == "Copy");
    BOOST_TEST(mgr.m_copy_as_menu->title().toStdString() == "Copy As");

    // All formats are present as actions, we just show/hide based on whether
    // there is a mol or a reaction present
    BOOST_TEST(mgr.m_copy_as_menu->actions().size() ==
               get_standard_export_formats().size() +
                   get_reaction_export_formats().size() +
                   2); // + (separator + image)

    // confirm copy as menu toggles based on reactions
    auto reaction_actions_visible = [&mgr](bool expect_reaction) {
        for (auto act : mgr.m_copy_as_menu->actions()) {
            bool expected = act->data().toBool() == expect_reaction;
            // copy as image (and its separator) are only visible for Copy All
            // As
            if (act->isSeparator() || act->text() == "Image") {
                expected =
                    mgr.m_copy_as_menu->title().toStdString() == "Copy All As";
            }
            if (act->isVisible() != expected) {
                return false;
            }
        }
        return true;
    };

    auto select_one_atom = [mol_model]() {
        auto atom = mol_model->getMol()->getAtomWithIdx(0);
        mol_model->select({atom}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    };

    BOOST_TEST(reaction_actions_visible(false));
    mol_model->clear();
    BOOST_TEST(reaction_actions_visible(false));
    sk.addFromString("CC>>CC");
    BOOST_TEST(reaction_actions_visible(true));
    select_one_atom();
    BOOST_TEST(reaction_actions_visible(false));
    mol_model->selectAll();
    BOOST_TEST(reaction_actions_visible(true));
    mol_model->clear();
    BOOST_TEST(reaction_actions_visible(false));
    sk.addFromString("CC");
    BOOST_TEST(reaction_actions_visible(false));

    // background context menu always "Copy All"
    mgr.setAlwaysCopyAll(true);

    BOOST_TEST(reaction_actions_visible(false));
    mol_model->clear();
    BOOST_TEST(reaction_actions_visible(false));
    sk.addFromString("CC>>CC");
    BOOST_TEST(reaction_actions_visible(true));
    select_one_atom();
    BOOST_TEST(reaction_actions_visible(true)); // ignores selection
    mol_model->selectAll();
    BOOST_TEST(reaction_actions_visible(true)); // ignores selection
    mol_model->clear();
    BOOST_TEST(reaction_actions_visible(false));
    sk.addFromString("CC");
    BOOST_TEST(reaction_actions_visible(false));
}

} // namespace sketcher
} // namespace schrodinger
