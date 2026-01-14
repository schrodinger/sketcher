#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/attachment_point_context_menu.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * SKETCH-2556: Verify that the attachment point context menu is created
 * correctly and has only the Delete action.
 */
BOOST_AUTO_TEST_CASE(test_attachment_point_menu_actions)
{
    AttachmentPointContextMenu menu;

    // The attachment point menu should have only one action: Delete
    auto actions = menu.actions();
    BOOST_TEST(actions.size() == 1);
    BOOST_TEST(actions[0]->text() == "Delete");
    BOOST_TEST(actions[0]->isEnabled() == true);
}

/**
 * SKETCH-2556: Verify that the menu can be populated with attachment point
 * atoms and bonds without errors, including bond-only selection.
 */
BOOST_AUTO_TEST_CASE(test_attachment_point_menu_set_context_items)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    AttachmentPointContextMenu menu;

    // Create a carbon atom with an attachment point
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(1);
    const auto* ap_bond = mol_model.getMol()->getBondWithIdx(0);

    BOOST_TEST(is_attachment_point(ap_atom));
    BOOST_TEST(is_attachment_point_bond(ap_bond));

    // Test 1: Set context to both atom and bond
    std::unordered_set<const RDKit::Atom*> atoms = {ap_atom};
    std::unordered_set<const RDKit::Bond*> bonds = {ap_bond};
    menu.setContextItems(atoms, bonds, {}, {});
    auto actions = menu.actions();
    BOOST_TEST(actions[0]->isEnabled() == true);

    // Test 2: Set context to only the bond (Scenario 2 from JIRA)
    menu.setContextItems({}, {ap_bond}, {}, {});
    actions = menu.actions();
    BOOST_TEST(actions[0]->isEnabled() == true);

    // Test 3: Set context to only the atom (Scenario 1 from JIRA)
    menu.setContextItems({ap_atom}, {}, {}, {});
    actions = menu.actions();
    BOOST_TEST(actions[0]->isEnabled() == true);
}

/**
 * SKETCH-2556: Verify that clicking Delete emits the deleteRequested signal
 * with attachment points.
 */
BOOST_AUTO_TEST_CASE(test_attachment_point_menu_delete_signal)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    AttachmentPointContextMenu menu;

    // Create a carbon atom with an attachment point
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(1);
    const auto* ap_bond = mol_model.getMol()->getBondWithIdx(0);

    std::unordered_set<const RDKit::Atom*> atoms = {ap_atom};
    std::unordered_set<const RDKit::Bond*> bonds = {ap_bond};
    menu.setContextItems(atoms, bonds, {}, {});

    // Connect to the deleteRequested signal and verify it's emitted
    bool signal_emitted = false;
    size_t emitted_atom_count = 0;
    size_t emitted_bond_count = 0;

    QObject::connect(
        &menu, &AttachmentPointContextMenu::deleteRequested,
        [&](const std::unordered_set<const RDKit::Atom*>& emitted_atoms,
            const std::unordered_set<const RDKit::Bond*>& emitted_bonds) {
            signal_emitted = true;
            emitted_atom_count = emitted_atoms.size();
            emitted_bond_count = emitted_bonds.size();
        });

    // Trigger the Delete action
    auto actions = menu.actions();
    BOOST_REQUIRE(actions.size() == 1);
    actions[0]->trigger();

    // Verify the signal was emitted with the expected number of items
    BOOST_TEST(signal_emitted);
    BOOST_TEST(emitted_atom_count == 1);
    BOOST_TEST(emitted_bond_count == 1);
}

} // namespace sketcher
} // namespace schrodinger
