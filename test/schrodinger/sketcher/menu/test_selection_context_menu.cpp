#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestSelectionContextMenu : public SelectionContextMenu
{
  public:
    TestSelectionContextMenu(SketcherModel* model, MolModel* mol_model) :
        SelectionContextMenu(model, mol_model)
    {
    }
    using SelectionContextMenu::m_clean_up_region_action;
    using SelectionContextMenu::m_variable_bond_action;
    using SelectionContextMenu::updateActions;
};

/**
 * Verify that the menu is properly updated to match the state of the model.
 */
BOOST_AUTO_TEST_CASE(test_updateActions)
{
    SketcherModel model;
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    auto menu = TestSelectionContextMenu(&model, &mol_model);
    import_mol_text(&mol_model, "NN.CCC");

    std::unordered_set<const RDKit::Atom*> c_atoms;
    std::unordered_set<const RDKit::Atom*> n_atoms;
    for (auto atom : mol_model.getMol()->atoms()) {
        switch (atom->getAtomicNum()) {
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

    // Select nothing
    menu.setContextItems({}, {}, {}, {}, {});
    // Use updateActions() as a placeholder for showEvent()
    menu.updateActions();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());

    // Select 1 atom
    menu.setContextItems({*c_atoms.begin()}, {}, {}, {}, {});
    menu.updateActions();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());

    // Select all atoms on the same molecule
    menu.setContextItems(n_atoms, {}, {}, {}, {});
    menu.updateActions();
    BOOST_TEST(menu.m_variable_bond_action->isEnabled());

    menu.setContextItems(c_atoms, {}, {}, {}, {});
    menu.updateActions();
    BOOST_TEST(menu.m_variable_bond_action->isEnabled());

    // Select atoms from different molecules
    menu.setContextItems({*c_atoms.begin(), *n_atoms.begin()}, {}, {}, {}, {});
    menu.updateActions();
    BOOST_TEST(!menu.m_variable_bond_action->isEnabled());

    // clean up region action should be enabled if there is a contiguous
    // region of atoms and bonds
    menu.setContextItems(c_atoms,
                         {mol_model.getMol()->getBondWithIdx(1),
                          mol_model.getMol()->getBondWithIdx(2)},
                         {}, {}, {});
    menu.updateActions();
    BOOST_TEST(menu.m_clean_up_region_action->isEnabled());

    // clean up region action should be disabled if there is no contiguous
    // region of atoms and bonds
    menu.setContextItems(c_atoms, {mol_model.getMol()->getBondWithIdx(0)}, {},
                         {}, {});
    menu.updateActions();
    BOOST_TEST(!menu.m_clean_up_region_action->isEnabled());
}

} // namespace sketcher
} // namespace schrodinger
