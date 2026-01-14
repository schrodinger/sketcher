
#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/ROMol.h>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/bracket_subgroup_context_menu.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestBracketSubgroupContextMenu : public BracketSubgroupContextMenu
{
  public:
    TestBracketSubgroupContextMenu() : BracketSubgroupContextMenu(nullptr){};
    using BracketSubgroupContextMenu::m_modify_notation_action;
    using BracketSubgroupContextMenu::setContextItems;
    using BracketSubgroupContextMenu::updateActions;
};

/**
 * Verify that actions are enabled/disabled as expected.
 */
BOOST_AUTO_TEST_CASE(test_updateActions)
{
    std::string smiles = "CC |SgD:0:FOO:42::::,SgD:1:BAR:3.14::::|";
    auto mol = rdkit_extensions::to_rdkit(smiles);
    auto sgroups = RDKit::getSubstanceGroups(*mol);
    BOOST_REQUIRE(sgroups.size() == 2);

    // modify should only be enabled for a single sgroup
    TestBracketSubgroupContextMenu menu;
    menu.setContextItems({}, {}, {}, {}, {});
    menu.updateActions();
    BOOST_TEST(!menu.m_modify_notation_action->isEnabled());

    menu.setContextItems({}, {}, {}, {&sgroups[0]}, {});
    menu.updateActions();
    BOOST_TEST(menu.m_modify_notation_action->isEnabled());

    menu.setContextItems({}, {}, {}, {&sgroups[0], &sgroups[1]}, {});
    menu.updateActions();
    BOOST_TEST(!menu.m_modify_notation_action->isEnabled());
}

} // namespace sketcher
} // namespace schrodinger