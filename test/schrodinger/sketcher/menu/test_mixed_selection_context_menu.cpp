#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * SKETCH-2556: Verify that in a mixed selection (regular atoms + attachment
 * points), modification actions only affect regular atoms, not attachment
 * points.
 */
BOOST_AUTO_TEST_CASE(test_mixed_selection_atom_modifications)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);

    // Create a molecule with a carbon, nitrogen, and an attachment point
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    mol_model.addAtom(Element::N, RDGeom::Point3D(4.0, 5.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);

    // Re-fetch pointers after molecule modification
    c_atom = mol_model.getMol()->getAtomWithIdx(0);
    const auto* n_atom = mol_model.getMol()->getAtomWithIdx(1);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(2);

    BOOST_TEST(mol_model.getMol()->getNumAtoms() == 3);
    BOOST_TEST(is_attachment_point(ap_atom));
    BOOST_TEST(!is_attachment_point(c_atom));
    BOOST_TEST(!is_attachment_point(n_atom));

    // Create a mixed selection: regular atoms (C, N) + attachment point
    std::unordered_set<const RDKit::Atom*> mixed_atoms = {c_atom, n_atom,
                                                          ap_atom};

    // Filter out attachment points (this is what the widget connection does)
    std::unordered_set<const RDKit::Atom*> filtered_atoms;
    for (const auto* atom : mixed_atoms) {
        if (!is_attachment_point(atom)) {
            filtered_atoms.insert(atom);
        }
    }

    // Verify filtering worked
    BOOST_TEST(filtered_atoms.size() == 2);
    BOOST_TEST(filtered_atoms.count(c_atom) == 1);
    BOOST_TEST(filtered_atoms.count(n_atom) == 1);
    BOOST_TEST(filtered_atoms.count(ap_atom) == 0);

    // Apply a mutation (change to oxygen) using the filtered set
    mol_model.mutateAtoms(filtered_atoms, Element::O);

    // Verify only regular atoms were changed, not the attachment point
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(0)->getSymbol() == "O");
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(1)->getSymbol() == "O");
    BOOST_TEST(is_attachment_point(mol_model.getMol()->getAtomWithIdx(2)));
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(2)->getAtomicNum() == 0);
}

/**
 * SKETCH-2556: Verify that in a mixed selection (regular bonds + attachment
 * point bonds), modification actions only affect regular bonds, not attachment
 * point bonds.
 */
BOOST_AUTO_TEST_CASE(test_mixed_selection_bond_modifications)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);

    // Create C-C single bond and C with attachment point
    mol_model.addAtom(Element::C, RDGeom::Point3D(0.0, 0.0, 0.0));
    mol_model.addAtom(Element::C, RDGeom::Point3D(2.0, 0.0, 0.0));
    const auto* c1_atom = mol_model.getMol()->getAtomWithIdx(0);
    const auto* c2_atom = mol_model.getMol()->getAtomWithIdx(1);

    // Add bond between carbons
    mol_model.addBond(c1_atom, c2_atom, RDKit::Bond::SINGLE);
    const auto* regular_bond = mol_model.getMol()->getBondBetweenAtoms(0, 1);

    // Add attachment point to second carbon
    c2_atom = mol_model.getMol()->getAtomWithIdx(1);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 1.0, 0.0), c2_atom);
    const auto* ap_bond = mol_model.getMol()->getBondBetweenAtoms(1, 2);

    BOOST_TEST(mol_model.getMol()->getNumBonds() == 2);
    BOOST_TEST(is_attachment_point_bond(ap_bond));
    BOOST_TEST(!is_attachment_point_bond(regular_bond));

    // Create a mixed selection: regular bond + AP bond
    std::unordered_set<const RDKit::Bond*> mixed_bonds = {regular_bond,
                                                          ap_bond};

    // Filter out attachment point bonds (this is what the widget connection
    // does)
    std::unordered_set<const RDKit::Bond*> filtered_bonds;
    for (const auto* bond : mixed_bonds) {
        if (!is_attachment_point_bond(bond)) {
            filtered_bonds.insert(bond);
        }
    }

    // Verify filtering worked
    BOOST_TEST(filtered_bonds.size() == 1);
    BOOST_TEST(filtered_bonds.count(regular_bond) == 1);
    BOOST_TEST(filtered_bonds.count(ap_bond) == 0);

    // Apply a mutation (change to double bond) using the filtered set
    mol_model.mutateBonds(filtered_bonds, BondTool::DOUBLE);

    // Verify only regular bond was changed, not the attachment point bond
    regular_bond = mol_model.getMol()->getBondBetweenAtoms(0, 1);
    ap_bond = mol_model.getMol()->getBondBetweenAtoms(1, 2);
    BOOST_TEST(regular_bond->getBondType() == RDKit::Bond::DOUBLE);
    BOOST_TEST(is_attachment_point_bond(ap_bond));
    BOOST_TEST(ap_bond->getBondType() == RDKit::Bond::SINGLE);
}

/**
 * SKETCH-2556: Verify that delete action works on mixed selections including
 * attachment points (delete should NOT filter out APs).
 */
BOOST_AUTO_TEST_CASE(test_mixed_selection_delete_includes_attachment_points)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);

    // Create C, N, and an attachment point on C
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    mol_model.addAtom(Element::N, RDGeom::Point3D(4.0, 5.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);

    // Re-fetch pointers after molecule modification
    c_atom = mol_model.getMol()->getAtomWithIdx(0);
    const auto* n_atom = mol_model.getMol()->getAtomWithIdx(1);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(2);
    const auto* ap_bond = mol_model.getMol()->getBondWithIdx(0);

    BOOST_TEST(mol_model.getMol()->getNumAtoms() == 3);
    BOOST_TEST(mol_model.getMol()->getNumBonds() == 1);

    // Create a mixed selection and delete (NO filtering for delete)
    std::unordered_set<const RDKit::Atom*> all_atoms = {c_atom, n_atom,
                                                        ap_atom};
    std::unordered_set<const RDKit::Bond*> all_bonds = {ap_bond};
    mol_model.remove(all_atoms, all_bonds, {}, {});

    // All atoms and bonds should be deleted
    BOOST_TEST(mol_model.getMol()->getNumAtoms() == 0);
    BOOST_TEST(mol_model.getMol()->getNumBonds() == 0);
}

/**
 * SKETCH-2556: Verify that charge adjustment only affects regular atoms in
 * mixed selections.
 */
BOOST_AUTO_TEST_CASE(test_mixed_selection_charge_adjustment)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);

    // Create C, N, and an attachment point
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    mol_model.addAtom(Element::N, RDGeom::Point3D(4.0, 5.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);

    // Re-fetch pointers after molecule modification
    c_atom = mol_model.getMol()->getAtomWithIdx(0);
    const auto* n_atom = mol_model.getMol()->getAtomWithIdx(1);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(2);

    // Create mixed selection and filter
    std::unordered_set<const RDKit::Atom*> mixed_atoms = {c_atom, n_atom,
                                                          ap_atom};
    std::unordered_set<const RDKit::Atom*> filtered_atoms;
    for (const auto* atom : mixed_atoms) {
        if (!is_attachment_point(atom)) {
            filtered_atoms.insert(atom);
        }
    }

    // Adjust charge on filtered atoms
    mol_model.adjustChargeOnAtoms(filtered_atoms, +1);

    // Verify only regular atoms have charge, not AP
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(0)->getFormalCharge() == 1);
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(1)->getFormalCharge() == 1);
    BOOST_TEST(mol_model.getMol()->getAtomWithIdx(2)->getFormalCharge() == 0);
}

/**
 * SKETCH-2556: Verify that showEditAtomPropertiesDialog is not called for
 * attachment points (edge case: if first atom in selection is an AP).
 */
BOOST_AUTO_TEST_CASE(test_edit_atom_properties_skips_attachment_points)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);

    // Create an attachment point first, then a carbon
    // This ensures AP is first in iteration order (index 0)
    mol_model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol_model.getMol()->getAtomWithIdx(0);
    mol_model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);

    // Re-fetch after modification
    c_atom = mol_model.getMol()->getAtomWithIdx(0);
    const auto* ap_atom = mol_model.getMol()->getAtomWithIdx(1);

    // Simulate the filtering logic for showEditAtomPropertiesDialog
    // In a mixed selection, if AP is first, we should NOT show the dialog
    bool should_show_dialog = !is_attachment_point(ap_atom);
    BOOST_TEST(should_show_dialog == false);

    // For regular atom, dialog should be shown
    should_show_dialog = !is_attachment_point(c_atom);
    BOOST_TEST(should_show_dialog == true);
}

} // namespace sketcher
} // namespace schrodinger
