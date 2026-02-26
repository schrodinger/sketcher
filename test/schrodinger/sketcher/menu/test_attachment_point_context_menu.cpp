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
 * SKETCH-2556: Verify that clicking Delete emits the deleteRequested signal
 * with the correct attachment point atoms and bonds.
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

    // Scenario 1: atom + bond context (right-clicking the attachment point
    // atom)
    std::unordered_set<const RDKit::Atom*> emitted_atoms;
    std::unordered_set<const RDKit::Bond*> emitted_bonds;
    QObject::connect(&menu, &AttachmentPointContextMenu::deleteRequested,
                     [&](const std::unordered_set<const RDKit::Atom*>& atoms,
                         const std::unordered_set<const RDKit::Bond*>& bonds) {
                         emitted_atoms = atoms;
                         emitted_bonds = bonds;
                     });

    menu.setContextItems({ap_atom}, {ap_bond}, {}, {}, {});
    menu.actions()[0]->trigger();
    BOOST_TEST(emitted_atoms ==
               std::unordered_set<const RDKit::Atom*>{ap_atom});
    BOOST_TEST(emitted_bonds ==
               std::unordered_set<const RDKit::Bond*>{ap_bond});

    // Scenario 2: bond-only context (right-clicking the attachment point bond)
    emitted_atoms.clear();
    emitted_bonds.clear();
    menu.setContextItems({}, {ap_bond}, {}, {}, {});
    menu.actions()[0]->trigger();
    BOOST_TEST(emitted_atoms.empty());
    BOOST_TEST(emitted_bonds ==
               std::unordered_set<const RDKit::Bond*>{ap_bond});
}

} // namespace sketcher
} // namespace schrodinger
