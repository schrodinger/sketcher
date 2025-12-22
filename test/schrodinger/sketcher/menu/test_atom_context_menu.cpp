
#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestModifyAtomsMenu : public ModifyAtomsMenu
{
  public:
    TestModifyAtomsMenu(SketcherModel* model, MolModel* mol_model) :
        ModifyAtomsMenu(model, mol_model){};
    using ModifyAtomsMenu::updateActions;
};

/**
 * @return The number of actions on this menu that are disabled
 */
int get_disabled_action_count(TestModifyAtomsMenu& menu)
{
    menu.updateActions();
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
    SketcherModel model;
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    auto menu = TestModifyAtomsMenu(&model, &mol_model);

    // Element, query, and an R group
    auto mdl = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.885714 -0.771429 0.000000 0
M  V30 2 A -3.648535 -0.057143 0.000000 0
M  V30 3 R# -2.411356 -0.771429 0.000000 0 RGROUPS=(1 1)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
)";
    auto mol = rdkit_extensions::to_rdkit(mdl);
    auto c_atom = mol->getAtomWithIdx(0);
    auto a_atom = mol->getAtomWithIdx(1);
    auto r_atom = mol->getAtomWithIdx(2);
    std::unordered_set<const RDKit::Atom*> atoms = {c_atom, a_atom, r_atom};

    // Without anything selected, the actions should all be disabled
    menu.setContextItems({}, {}, {}, {}, {});
    int disabled_count = get_disabled_action_count(menu);
    BOOST_TEST(disabled_count == 6);

    // Try specifying individual atoms: most actions should be enabled
    for (auto atom : atoms) {
        int exp_disabled_count = 0;

        if (is_r_group(atom)) {
            // The R Group atom type should have all actions disabled
            exp_disabled_count = 5;
        } else if (atom->hasQuery()) {
            // The Any Atom Wildcard type should only have setCharge enabled
            exp_disabled_count = 3;
        } else if (atom->getAtomicNum() > 0) {
            // Only real atoms should have all options except Add Brackets
            exp_disabled_count = 1;
        }

        menu.setContextItems({atom}, {}, {}, {}, {});
        disabled_count = get_disabled_action_count(menu);
        BOOST_TEST(disabled_count == exp_disabled_count);
    }

    // Test after altering the state of the element atom
    menu.setContextItems({c_atom}, {}, {}, {}, {});
    for (unsigned int unpaired_e_count = MIN_UNPAIRED_E;
         unpaired_e_count <= MAX_UNPAIRED_E; ++unpaired_e_count) {
        c_atom->setNumRadicalElectrons(unpaired_e_count);
        bool unpaired_at_extreme =
            unpaired_e_count == 0 || unpaired_e_count == 4;

        for (int charge = -ATOM_CHARGE_LIMIT; charge <= ATOM_CHARGE_LIMIT;
             ++charge) {

            bool charge_at_extreme = std::abs(charge) == ATOM_CHARGE_LIMIT;

            c_atom->setFormalCharge(charge);
            int exp_disabled_count =
                1 + unpaired_at_extreme + charge_at_extreme;

            // Assign all objects in the scene to the context menu
            menu.setContextItems(atoms, {}, {}, {}, {});
            disabled_count = get_disabled_action_count(menu);
            BOOST_TEST(disabled_count == exp_disabled_count);

            menu.setContextItems({c_atom}, {}, {}, {}, {});
            --exp_disabled_count;
            disabled_count = get_disabled_action_count(menu);
            BOOST_TEST(disabled_count == exp_disabled_count);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger