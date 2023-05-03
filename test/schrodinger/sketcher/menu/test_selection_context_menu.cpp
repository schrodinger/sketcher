#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "../test_sketcherScene.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

class TestSelectionContextMenu : public SelectionContextMenu
{
  public:
    TestSelectionContextMenu(SketcherModel* model) : SelectionContextMenu(model)
    {
    }
    using SelectionContextMenu::m_variable_bond_action;
    using SelectionContextMenu::updateActionsEnabled;
};

/**
 * Verify that the menu is properly updated to match the state of the model.
 */
BOOST_AUTO_TEST_CASE(test_updateActionsEnabled)
{
    testSketcherScene scene;
    scene.importText("NN.CCC");
    auto model = scene.getModel();
    TestSelectionContextMenu menu(model);
    scene.connectContextMenu(menu);

    std::unordered_set<sketcherAtom*> c_atoms;
    std::unordered_set<sketcherAtom*> n_atoms;
    for (auto atom : scene.quickGetAtoms()) {
        switch (atom->getAtomType()) {
            case static_cast<int>(Element::C): {
                c_atoms.insert(atom);
                break;
            }
            case static_cast<int>(Element::N): {
                n_atoms.insert(atom);
                break;
            }
        }
    }

    // Variable bond will be disabled if fewer than 2 atoms are selected
    std::unordered_set<sketcherAtom*> selection;

    // Select nothing
    scene.setSelection(selection);
    scene.clearContextMenuObjects();
    menu.updateActionsEnabled();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());

    // Select 1 atom
    selection.insert(*c_atoms.begin());
    scene.setSelection(selection);
    scene.clearContextMenuObjects();
    for (const auto& obj : scene.getSelectedObjects()) {
        scene.m_context_menu_objects.insert(obj);
    }
    menu.updateActionsEnabled();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());

    // Select all atoms on the same molecule
    scene.setSelection(n_atoms);
    scene.clearContextMenuObjects();
    for (const auto& obj : scene.getSelectedObjects()) {
        scene.m_context_menu_objects.insert(obj);
    }
    menu.updateActionsEnabled();
    BOOST_TEST(menu.m_variable_bond_action->isEnabled());

    scene.setSelection(c_atoms);
    scene.clearContextMenuObjects();
    for (const auto& obj : scene.getSelectedObjects()) {
        scene.m_context_menu_objects.insert(obj);
    }
    menu.updateActionsEnabled();
    BOOST_TEST(menu.m_variable_bond_action->isEnabled());

    // Select atoms from different molecules
    selection.clear();
    selection.insert(*c_atoms.begin());
    selection.insert(*n_atoms.begin());
    scene.setSelection(selection);
    scene.clearContextMenuObjects();
    for (const auto& obj : scene.getSelectedObjects()) {
        scene.m_context_menu_objects.insert(obj);
    }
    menu.updateActionsEnabled();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());
}

} // namespace sketcher
} // namespace schrodinger
