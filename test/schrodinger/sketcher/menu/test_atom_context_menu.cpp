#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "../test_sketcherScene.h"
#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * @return The number of actions on this menu that are disabled
 */
int get_disabled_action_count(ModifyAtomsMenu& menu)
{
    menu.updateActionsEnabled();
    int disabled_count = 0;
    for (auto action : menu.actions()) {
        if (!action->isEnabled()) {
            ++disabled_count;
        }
    }
    return disabled_count;
}

/**
 * Verify that actions are enabled/disabled as expected.
 */
BOOST_AUTO_TEST_CASE(test_actions_enabled)
{

    testSketcherScene scene;
    auto menu = ModifyAtomsMenu(scene.getModel());

    scene.importText("CCC");
    std::vector<sketcherAtom*> atoms;
    scene.getAtoms(atoms);
    auto c_atom = atoms[0];
    auto a_atom = atoms[1];
    auto r_atom = atoms[2];

    // Change one of the atoms to a query
    a_atom->setAtomType(AtomTypes::A_QUERY_KEY);

    // Replace another atom with an R group
    r_atom->setAtomType(R_GROUP_KEY);

    // Without anything selected, the actions should all be disabled
    scene.clearContextMenuObjects();
    int disabled_count = get_disabled_action_count(menu);
    BOOST_TEST(disabled_count == 6);

    // Try specifying individual atoms: most actions should be enabled
    for (auto atom : scene.quickGetAtoms()) {
        int exp_disabled_count = 0;

        if (is_atomic_number(atom->getAtomType())) {
            // Only real atoms should have all options except Add Brackets
            exp_disabled_count = 1;
        } else if (atom->getAtomType() == A_QUERY_KEY) {
            // The Any Atom Wildcard type should only have setCharge enabled
            exp_disabled_count = 3;
        } else if (atom->getAtomType() == R_GROUP_KEY) {
            // The R Group atom type should have all actions disabled
            exp_disabled_count = 5;
        }

        scene.m_context_menu_objects.clear();
        scene.m_context_menu_objects.insert(atom);
        disabled_count = get_disabled_action_count(menu);
        BOOST_TEST(disabled_count == exp_disabled_count);
    }

    // Test after altering the state of the element atom
    scene.clearContextMenuObjects();
    scene.m_context_menu_objects.insert(c_atom);
    for (unsigned int unpaired_e_count = MIN_UNPAIRED_E;
         unpaired_e_count <= MAX_UNPAIRED_E; ++unpaired_e_count) {
        c_atom->setUnpairedElectronsN(unpaired_e_count);
        bool unpaired_at_extreme =
            unpaired_e_count == 0 || unpaired_e_count == 4;

        for (int charge = -ATOM_CHARGE_LIMIT; charge <= ATOM_CHARGE_LIMIT;
             ++charge) {

            bool charge_at_extreme = std::abs(charge) == ATOM_CHARGE_LIMIT;

            c_atom->setCharge(charge);
            int exp_disabled_count =
                1 + unpaired_at_extreme + charge_at_extreme;

            // Assign all objects in the scene to the context menu
            for (const auto& obj : scene.getObjects()) {
                scene.m_context_menu_objects.insert(obj);
            }

            disabled_count = get_disabled_action_count(menu);
            BOOST_TEST(disabled_count == exp_disabled_count);
            scene.clearContextMenuObjects();
            scene.m_context_menu_objects.insert(c_atom);
            --exp_disabled_count;
            disabled_count = get_disabled_action_count(menu);
            BOOST_TEST(disabled_count == exp_disabled_count);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger