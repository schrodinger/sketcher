#define BOOST_TEST_MODULE Test_Sketcher

#include "schrodinger/sketcher/molviewer/mol_model.h"

#include <unordered_set>

#include <QUndoStack>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "../test_common.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);
// Boost doesn't know how to print QStrings
BOOST_TEST_DONT_PRINT_LOG_VALUE(QString);
BOOST_TEST_DONT_PRINT_LOG_VALUE(RDGeom::Point3D);

namespace schrodinger
{
namespace sketcher
{

typedef std::unordered_set<const RDKit::Atom*> atom_set;
typedef std::unordered_set<const RDKit::Bond*> bond_set;

class TestMolModel : public MolModel
{
  public:
    TestMolModel(QUndoStack* const undo_stack) : MolModel(undo_stack)
    {
    }
    using MolModel::getAtomFromTag;
    using MolModel::getBondFromTag;
    using MolModel::getTagForAtom;
    using MolModel::getTagForBond;
};

void check_coords(const RDKit::ROMol* mol, unsigned int atom_index,
                  double exp_x, double exp_y)
{
    const RDGeom::Point3D& point = mol->getConformer().getAtomPos(atom_index);
    double difference = (point - RDGeom::Point3D(exp_x, exp_y, 0.0)).lengthSq();
    BOOST_TEST(difference < 0.001);
}

BOOST_AUTO_TEST_CASE(test_addAtom)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 0);
    model.addAtom("C", RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    check_coords(mol, 0, 1.0, 2.0);
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 2);
    auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    // make sure c_atom is still valid
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == c_atom);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);

    // undo again to make sure that redoing hasn't modified the structures
    // stored in the undo lambda
    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == c_atom);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
}

BOOST_AUTO_TEST_CASE(test_removeAtom)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom("C", RDGeom::Point3D(1.0, 2.0, 0.0));
    auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    auto* n_atom = mol->getAtomWithIdx(1);

    model.removeAtom(c_atom);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == n_atom);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    model.removeAtom(n_atom);
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    n_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(model.getAtomFromTag(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == n_atom);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
}

BOOST_AUTO_TEST_CASE(test_addBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom("C", RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    const auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(mol->getNumBonds() == 0);

    model.addBond(c_atom, n_atom, RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(0) == bond);
    BOOST_TEST(model.getTagForBond(bond) == 0);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(0) == bond);
    BOOST_TEST(model.getTagForBond(bond) == 0);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
}

BOOST_AUTO_TEST_CASE(test_removeBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom("C", RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    const auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(mol->getNumBonds() == 0);

    model.addBond(c_atom, n_atom, RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* bond = mol->getBondWithIdx(0);

    model.removeBond(bond);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(0) == bond);
    BOOST_TEST(model.getTagForBond(bond) == 0);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 0);
}

BOOST_AUTO_TEST_CASE(test_addMol)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);

    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCC"));
    model.addMol(*mol_to_add);
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
    const RDKit::Atom* atom;
    for (int i = 1; i < 5; ++i) {
        atom = mol->getAtomWithIdx(i);
        BOOST_TEST(atom->getSymbol() == "C");
        BOOST_TEST(model.getAtomFromTag(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(i) == bond);
        BOOST_TEST(model.getTagForBond(bond) == i);
    }

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    BOOST_TEST(mol->getAtomWithIdx(0)->getSymbol() == "N");

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
    for (int i = 1; i < 5; ++i) {
        atom = mol->getAtomWithIdx(i);
        BOOST_TEST(atom->getSymbol() == "C");
        BOOST_TEST(model.getAtomFromTag(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(i) == bond);
        BOOST_TEST(model.getTagForBond(bond) == i);
    }
}

BOOST_AUTO_TEST_CASE(test_clear)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    model.clear();
    BOOST_TEST(!undo_stack.count());

    const RDKit::ROMol* mol = model.getMol();
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCC"));
    model.addMol(*mol_to_add);
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);

    model.clear();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
    const RDKit::Atom* atom;
    for (int i = 1; i < 5; ++i) {
        atom = mol->getAtomWithIdx(i);
        BOOST_TEST(atom->getSymbol() == "C");
        BOOST_TEST(model.getAtomFromTag(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(i) == bond);
        BOOST_TEST(model.getTagForBond(bond) == i);
    }

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
}

BOOST_AUTO_TEST_CASE(test_selection)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());

    // calls that don't change anything shouldn't add a command to the undo
    // stack
    model.clearSelection();
    BOOST_TEST(!undo_stack.count());
    model.select({}, {}, SelectMode::SELECT);
    BOOST_TEST(!undo_stack.count());
    model.select({}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(!undo_stack.count());

    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCCC"));
    model.addMol(*mol_to_add);
    BOOST_TEST(undo_stack.count() == 1);

    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());

    // calls that don't change anything shouldn't add a command to the undo
    // stack
    model.clearSelection();
    BOOST_TEST(undo_stack.count() == 1);
    model.select({}, {}, SelectMode::SELECT);
    BOOST_TEST(undo_stack.count() == 1);
    model.select({}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(undo_stack.count() == 1);

    auto* atom1 = model.getAtomFromTag(1);
    auto* atom2 = model.getAtomFromTag(2);
    auto* bond1 = model.getBondFromTag(1);
    model.select({atom1, atom2}, {bond1}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    model.select({atom1}, {bond1}, SelectMode::DESELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    model.clearSelection();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());

    model.select({atom2}, {bond1}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    // make sure that selecting atoms/bonds that are already selected can be
    // undone correctly
    auto* bond2 = model.getBondFromTag(2);
    model.select({atom1, atom2}, {bond1, bond2}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));

    // make sure that deselecting atoms/bonds that are already deselected can be
    // undone correctly
    auto* atom3 = model.getAtomFromTag(3);
    auto* bond3 = model.getBondFromTag(3);
    model.select({atom2, atom3}, {bond2, bond3}, SelectMode::DESELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    // toggle the selection
    model.select({atom1}, {}, SelectMode::TOGGLE);
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    model.select({atom1}, {}, SelectMode::TOGGLE);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    // select-only
    model.select({atom2}, {bond3}, SelectMode::SELECT_ONLY);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));

    // select-only with no atoms or bonds specified should clear the selection
    model.select({}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));
}

BOOST_AUTO_TEST_CASE(test_select_all_and_invert)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCCC"));
    model.addMol(*mol_to_add);

    model.selectAll();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);

    auto* atom1 = model.getAtomFromTag(1);
    auto* atom2 = model.getAtomFromTag(2);
    auto* bond1 = model.getBondFromTag(1);
    model.select({atom1, atom2}, {bond1}, SelectMode::SELECT_ONLY);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    model.selectAll();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    auto* atom3 = model.getAtomFromTag(3);
    auto* atom4 = model.getAtomFromTag(4);
    auto* atom0 = model.getAtomFromTag(0);
    auto* bond2 = model.getBondFromTag(2);
    auto* bond3 = model.getBondFromTag(3);
    auto* bond0 = model.getBondFromTag(0);
    model.invertSelection();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom3, atom4, atom0}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond2, bond3, bond0}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom3, atom4, atom0}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond2, bond3, bond0}));
}

/**
 * Ensure that removing a selected atom or bond automatically deselects the
 * removed element.  Ensure that undoing the removal restores the selection.
 */
BOOST_AUTO_TEST_CASE(test_removing_selected)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCCC"));
    model.addMol(*mol_to_add);
    auto* atom1 = model.getAtomFromTag(1);
    auto* atom2 = model.getAtomFromTag(2);
    auto* bond1 = model.getBondFromTag(1);
    model.select({atom1, atom2}, {bond1}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    model.removeBond(bond1);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());

    undo_stack.undo();
    atom1 = model.getAtomFromTag(1);
    atom2 = model.getAtomFromTag(2);
    bond1 = model.getBondFromTag(1);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    // removing an atom removes all of its bonds, so this removeAtom call will
    // implicitly remove and deselect bond1
    model.removeAtom(model.getAtomFromTag(1));
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());

    undo_stack.undo();
    atom1 = model.getAtomFromTag(1);
    atom2 = model.getAtomFromTag(2);
    bond1 = model.getBondFromTag(1);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    model.clear();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());

    undo_stack.undo();
    atom1 = model.getAtomFromTag(1);
    atom2 = model.getAtomFromTag(2);
    bond1 = model.getBondFromTag(1);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    // removing non-selected atoms and bonds shouldn't have any effect on the
    // selection
    model.removeAtom(model.getAtomFromTag(3));
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));

    undo_stack.undo();
    model.removeBond(model.getBondFromTag(2));
    atom1 = model.getAtomFromTag(1);
    atom2 = model.getAtomFromTag(2);
    bond1 = model.getBondFromTag(1);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
}

} // namespace sketcher
} // namespace schrodinger
