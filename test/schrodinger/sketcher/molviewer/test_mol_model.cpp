#define BOOST_TEST_MODULE Test_Sketcher

#include "schrodinger/sketcher/molviewer/mol_model.h"

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

class TestMolModel : public MolModel
{
  public:
    TestMolModel(QUndoStack* const undo_stack) : MolModel(undo_stack)
    {
    }
    using MolModel::getTagForAtom;
    using MolModel::getTagForBond;

    // RDKit is missing const versions of bookmark getters, so we use
    // const_cast
    const RDKit::Atom* getAtom(int atom_tag) const
    {
        return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueAtomWithBookmark(
            atom_tag);
    }

    const RDKit::Bond* getBond(int bond_tag) const
    {
        return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueBondWithBookmark(
            bond_tag);
    }
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
    BOOST_TEST(model.getAtom(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    check_coords(mol, 0, 1.0, 2.0);
    model.addAtom("N", RDGeom::Point3D(3.0, 4.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 2);
    auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtom(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    // make sure c_atom is still valid
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtom(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == c_atom);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtom(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtom(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtom(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtom(0) == c_atom);
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
    BOOST_TEST(model.getAtom(0) == c_atom);
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
    BOOST_TEST(model.getAtom(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    model.removeAtom(n_atom);
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    n_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtom(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(model.getAtom(0) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == 0);
    BOOST_TEST(model.getAtom(1) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == 1);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == n_atom);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtom(1) == n_atom);
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
    BOOST_TEST(model.getBond(0) == bond);
    BOOST_TEST(model.getTagForBond(bond) == 0);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBond(0) == bond);
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
    BOOST_TEST(model.getBond(0) == bond);
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
        BOOST_TEST(model.getAtom(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBond(i) == bond);
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
        BOOST_TEST(model.getAtom(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBond(i) == bond);
        BOOST_TEST(model.getTagForBond(bond) == i);
    }
}

BOOST_AUTO_TEST_CASE(test_clear)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
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
        BOOST_TEST(model.getAtom(i) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == i);
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBond(i) == bond);
        BOOST_TEST(model.getTagForBond(bond) == i);
    }

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
}

} // namespace sketcher
} // namespace schrodinger
