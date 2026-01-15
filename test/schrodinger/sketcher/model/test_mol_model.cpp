/**
 * Tests for the MolModel class.
 *
 * When writing tests, note that *any* change to the MolModel's molecule will
 * invalidate *all* Atom and Bond pointers that reference the MolModel's
 * molecule.  For example, flipping a bond will invalidate any pointers to the
 * flipped bond, any pointers to any other bond in the molecule, and any
 * Atom* pointers to the molecule.
 */

#define BOOST_TEST_MODULE Test_Sketcher

#include <optional>
#include <unordered_set>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>
#include <QSignalSpy>

#include <boost/algorithm/string.hpp>
#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/MolTransforms/MolTransforms.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/QueryBond.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <QUndoStack>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"
#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/fragment.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/monomer_connectors.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/s_group_constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(RDGeom::Point3D);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::NonMolecularObject);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::NonMolecularType);
BOOST_TEST_DONT_PRINT_LOG_VALUE(RDKit::Bond::BondDir);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::rdkit_extensions::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::rdkit_extensions::EnhancedStereo>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::sketcher::EnhancedStereo>)

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

using r_group_num_t = std::optional<unsigned int>;
static auto no_r_group_num = std::nullopt;

BOOST_TEST_DONT_PRINT_LOG_VALUE(r_group_num_t);
BOOST_TEST_DONT_PRINT_LOG_VALUE(decltype(no_r_group_num));

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

typedef std::unordered_set<const RDKit::Atom*> atom_set;
typedef std::unordered_set<const RDKit::Bond*> bond_set;
typedef std::unordered_set<const NonMolecularObject*> nmo_set;

class TestMolModel : public MolModel
{
  public:
    TestMolModel(QUndoStack* const undo_stack) : MolModel(undo_stack)
    {
    }
    using MolModel::flipBondStereo;
    using MolModel::getAllUnselectedTags;
    using MolModel::getAtomFromTag;
    using MolModel::getBondFromTag;
    using MolModel::getTagForAtom;
    using MolModel::getTagForBond;
    using MolModel::m_arrow;
    using MolModel::m_pluses;
    using MolModel::m_tag_to_non_molecular_object;
};

const std::string STRUC_WITH_RINGS = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 28 32 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -0.180969 -4.827934 0.000000 0
M  V30 2 C 0.684718 -4.327610 0.000000 0
M  V30 3 N 0.684341 -3.327343 0.000000 0
M  V30 4 C 1.549928 -2.827192 0.000000 0
M  V30 5 N 1.549625 -1.827197 0.000000 0
M  V30 6 C 2.416364 -3.326887 0.000000 0
M  V30 7 N 3.282051 -2.826562 0.000000 0
M  V30 8 C 4.148313 -3.326158 0.000000 0
M  V30 9 C 5.014000 -2.825833 0.000000 0
M  V30 10 N 5.880263 -3.325429 0.000000 0
M  V30 11 C 6.706370 -2.761842 0.000000 0
M  V30 12 C 6.483613 -1.787165 0.000000 0
M  V30 13 C 7.216485 -1.106732 0.000000 0
M  V30 14 C 8.172187 -1.401248 0.000000 0
M  V30 15 C 8.394845 -2.376098 0.000000 0
M  V30 16 C 7.661974 -3.056531 0.000000 0
M  V30 17 C 8.027705 -3.987299 0.000000 0
M  V30 18 C 7.527762 -4.853363 0.000000 0
M  V30 19 C 6.539072 -5.002846 0.000000 0
M  V30 20 C 6.316569 -5.977700 0.000000 0
M  V30 21 C 5.361167 -6.272597 0.000000 0
M  V30 22 C 4.628094 -5.592542 0.000000 0
M  V30 23 C 4.850424 -4.617588 0.000000 0
M  V30 24 C 5.805826 -4.322690 0.000000 0
M  V30 25 C 4.148443 -4.326053 0.000000 0
M  V30 26 N 3.282756 -4.826377 0.000000 0
M  V30 27 C 2.416667 -4.326881 0.000000 0
M  V30 28 N 1.550807 -4.827106 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 2 4 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 2 11 12
M  V30 12 1 12 13
M  V30 13 2 13 14
M  V30 14 1 14 15
M  V30 15 2 15 16
M  V30 16 1 11 16
M  V30 17 1 16 17
M  V30 18 2 17 18
M  V30 19 1 18 19
M  V30 20 2 19 20
M  V30 21 1 20 21
M  V30 22 2 21 22
M  V30 23 1 22 23
M  V30 24 2 23 24
M  V30 25 1 10 24
M  V30 26 1 19 24
M  V30 27 1 8 25
M  V30 28 2 25 26
M  V30 29 1 26 27
M  V30 30 1 6 27
M  V30 31 2 27 28
M  V30 32 1 2 28
M  V30 END BOND
M  V30 END CTAB
M  END
)";

/**
 * A molecule with a variable attachment bond
 */
const std::string VAR_ATTACH_MOL = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.065751 0.028571 0.000000 0
M  V30 2 C 1.302930 0.742857 0.000000 0
M  V30 3 C 1.302930 2.171429 0.000000 0
M  V30 4 C 0.065751 2.885714 0.000000 0
M  V30 5 C -1.171429 2.171429 0.000000 0
M  V30 6 C -1.171429 0.742857 0.000000 0
M  V30 7 C 1.398626 3.765751 0.000000 0
M  V30 8 * 0.327197 1.909982 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 1 8 7 ENDPTS=(3 2 3 4) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$

)";

const std::string AMIDE_FRAG_SMILES = "*CC(N)=O |$_AP1;;;;$|";
const std::string CYCLOHEXANE_MOL_SMILES = "C1CCCCC1";
const std::string CYCLOHEXANE_FRAG_SMILES = "*C1CCCCC1 |$_AP1;;;;;;$|";
const std::string BENZENE_MOL_SMILES = "C1=CC=CC=C1";
const std::string BENZENE_FRAG_SMILES = "*C1=CC=CC=C1 |$_AP1;;;;;;$|";
const std::string METHYLCYCLOHEXANE_MOL_SMILES = "CC1CCCCC1";
const std::string TWO_RINGS_SMILES = "C1=CC=C2C=CC=CC2=C1";

void check_coords(const RDGeom::Point3D& point, const double exp_x,
                  const double exp_y)
{
    double difference = (point - RDGeom::Point3D(exp_x, exp_y, 0.0)).lengthSq();
    BOOST_TEST(difference < 0.001);
}

void check_coords(const RDKit::ROMol* mol, unsigned int atom_index,
                  const double exp_x, const double exp_y)
{
    const RDGeom::Point3D& point = mol->getConformer().getAtomPos(atom_index);
    check_coords(point, exp_x, exp_y);
}

BOOST_AUTO_TEST_CASE(test_addAtom)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 0);
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));
    check_coords(mol, 0, 1.0, 2.0);
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));
    // make sure c_atom is still valid
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == c_atom);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 1, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));

    // undo again to make sure that redoing hasn't modified the structures
    // stored in the undo lambda
    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getAtomWithIdx(0) == c_atom);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
}

BOOST_AUTO_TEST_CASE(test_addAtomWithBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom(Element::N, RDGeom::Point3D(1.0, 2.0, 0.0), c_atom,
                  RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::NONE);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    c_atom = mol->getAtomWithIdx(0);
    const auto* n_atom = mol->getAtomWithIdx(1);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(BondTag(0)) == bond);
    BOOST_TEST(model.getTagForBond(bond) == BondTag(0));
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(BondTag(0)) == bond);
    BOOST_TEST(model.getTagForBond(bond) == BondTag(0));
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);
}

BOOST_AUTO_TEST_CASE(test_addAtomChain)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtomChain(Element::C,
                       {RDGeom::Point3D(1.0, 0.0, 0.0),
                        RDGeom::Point3D(2.0, 0.0, 0.0),
                        RDGeom::Point3D(3.0, 0.0, 0.0)},
                       nullptr, RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 2);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 2);

    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtomChain(
        Element::N,
        {RDGeom::Point3D(10.0, 1.0, 0.0), RDGeom::Point3D(1.0, 2.0, 0.0)},
        c_atom, RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::NONE);
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 4);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 2);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 4);
}

BOOST_AUTO_TEST_CASE(test_addAtom_query)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    auto atom_query = std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>(
        RDKit::makeAAtomQuery());
    model.addAtom(atom_query, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    auto* query_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(query_atom->hasQuery());
    BOOST_TEST(query_atom->getQueryType() == "A");

    auto bond_query = std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(
        RDKit::makeBondNullQuery());
    model.addAtom(Element::C, RDGeom::Point3D(3.0, 4.0, 0.0), query_atom,
                  bond_query);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    auto* c_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(!c_atom->hasQuery());
    BOOST_TEST(c_atom->getSymbol() == "C");
    auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->hasQuery());
    BOOST_TEST(bond->getQuery()->getTypeLabel() == bond_query->getTypeLabel());

    model.addAtom(atom_query, RDGeom::Point3D(5.0, 6.0, 0.0), c_atom,
                  bond_query);
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 2);
    auto* query_atom2 = mol->getAtomWithIdx(2);
    BOOST_TEST(query_atom2->hasQuery());
    BOOST_TEST(query_atom2->getQueryType() == "A");
    auto* bond2 = mol->getBondWithIdx(1);
    BOOST_TEST(bond2->hasQuery());
    BOOST_TEST(bond2->getQuery()->getTypeLabel() == bond_query->getTypeLabel());

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 2);
}

using schrodinger::rdkit_extensions::get_r_group_number;

BOOST_AUTO_TEST_CASE(test_addAtom_r_group)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    BOOST_TEST(get_all_r_group_numbers(mol).empty());
    BOOST_TEST(get_next_r_group_numbers(mol, 1) ==
               std::vector<unsigned int>({1}));
    BOOST_TEST(get_next_r_group_numbers(mol, 2) ==
               std::vector<unsigned int>({1, 2}));
    BOOST_TEST(get_next_r_group_numbers(mol, 3) ==
               std::vector<unsigned int>({1, 2, 3}));

    model.addRGroup(3, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    const auto* r3_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(r3_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r3_atom) == 3);
    BOOST_TEST(is_r_group(r3_atom));
    BOOST_TEST(!is_attachment_point(r3_atom));
    BOOST_TEST(get_all_r_group_numbers(mol) == std::vector<unsigned int>({3}));
    BOOST_TEST(get_next_r_group_numbers(mol, 1) ==
               std::vector<unsigned int>({1}));
    BOOST_TEST(get_next_r_group_numbers(mol, 2) ==
               std::vector<unsigned int>({1, 2}));
    BOOST_TEST(get_next_r_group_numbers(mol, 3) ==
               std::vector<unsigned int>({1, 2, 4}));

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    r3_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(r3_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r3_atom) == 3);
    BOOST_TEST(is_r_group(r3_atom));
    BOOST_TEST(!is_attachment_point(r3_atom));

    model.addRGroup(1, RDGeom::Point3D(3.0, 4.0, 0.0), r3_atom);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* r1_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(r1_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r1_atom) == 1);
    BOOST_TEST(is_r_group(r1_atom));
    BOOST_TEST(!is_attachment_point(r1_atom));
    BOOST_TEST(get_all_r_group_numbers(mol) ==
               std::vector<unsigned int>({1, 3}));
    BOOST_TEST(get_next_r_group_numbers(mol, 1) ==
               std::vector<unsigned int>({2}));
    BOOST_TEST(get_next_r_group_numbers(mol, 2) ==
               std::vector<unsigned int>({2, 4}));
    BOOST_TEST(get_next_r_group_numbers(mol, 3) ==
               std::vector<unsigned int>({2, 4, 5}));

    model.addRGroupChain({2, 4}, {RDGeom::Point3D(5.0, 6.0, 0.0),
                                  RDGeom::Point3D(7.0, 8.0, 0.0)});
    BOOST_TEST(mol->getNumAtoms() == 4);
    BOOST_TEST(mol->getNumBonds() == 2);
    const auto* r2_atom = mol->getAtomWithIdx(2);
    BOOST_TEST(r2_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r2_atom) == 2);
    const auto* r4_atom = mol->getAtomWithIdx(3);
    BOOST_TEST(r4_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r4_atom) == 4);
    BOOST_TEST(get_all_r_group_numbers(mol) ==
               std::vector<unsigned int>({1, 2, 3, 4}));
    BOOST_TEST(get_next_r_group_numbers(mol, 1) ==
               std::vector<unsigned int>({5}));
    BOOST_TEST(get_next_r_group_numbers(mol, 2) ==
               std::vector<unsigned int>({5, 6}));
    BOOST_TEST(get_next_r_group_numbers(mol, 3) ==
               std::vector<unsigned int>({5, 6, 7}));

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);

    undo_stack.redo();
    r2_atom = mol->getAtomWithIdx(2);
    BOOST_TEST(r2_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r2_atom) == 2);
    r4_atom = mol->getAtomWithIdx(3);
    BOOST_TEST(r4_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(r4_atom) == 4);

    // Confirm that the exported SMILES does not contain stray atom properties
    auto expected_smiles = "**.** |$_R3;_R1;_R2;_R4$|";
    auto export_smiles = rdkit_extensions::to_string(*model.getMolForExport(),
                                                     Format::EXTENDED_SMILES);
    BOOST_TEST(export_smiles == expected_smiles);

    // Confirm that imported rgroups roundtrip similarly
    model.clear();
    add_text_to_mol_model(model, expected_smiles);
    auto export_mol = model.getMolForExport();
    // Though these properties are present on rgroups....
    for (const auto* atom : export_mol->atoms()) {
        BOOST_TEST(atom->hasProp(RDKit::common_properties::atomLabel));
        BOOST_TEST(atom->hasProp(RDKit::common_properties::_MolFileRLabel));
    }
    // ...confirm the exported SMILES lack those properties in the extension
    export_smiles =
        rdkit_extensions::to_string(*export_mol, Format::EXTENDED_SMILES);
    BOOST_TEST(export_smiles == expected_smiles);
}

BOOST_AUTO_TEST_CASE(test_addAtom_attachment_point)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();

    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    const auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(!is_r_group(c_atom));
    BOOST_TEST(!is_attachment_point(c_atom));
    BOOST_TEST(get_r_group_number(c_atom) == no_r_group_num);
    BOOST_TEST(get_next_attachment_point_number(mol) == 1);

    model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* attachment_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(attachment_atom->getAtomicNum() == 0);
    BOOST_TEST(is_attachment_point(attachment_atom));
    BOOST_TEST(get_attachment_point_number(attachment_atom) == 1);
    BOOST_TEST(!is_r_group(attachment_atom));
    BOOST_TEST(get_r_group_number(attachment_atom) == no_r_group_num);
    BOOST_TEST(get_all_r_group_numbers(mol).empty());
    BOOST_TEST(get_next_r_group_numbers(mol, 1) ==
               std::vector<unsigned int>({1}));
    BOOST_TEST(get_next_attachment_point_number(mol) == 2);

    c_atom = mol->getAtomWithIdx(0);
    model.addAttachmentPoint(RDGeom::Point3D(3.0, 4.0, 0.0), c_atom);
    const auto* attachment_atom2 = mol->getAtomWithIdx(2);
    BOOST_TEST(attachment_atom2->getAtomicNum() == 0);
    BOOST_TEST(is_attachment_point(attachment_atom2));
    BOOST_TEST(get_attachment_point_number(attachment_atom2) == 2);
    BOOST_TEST(get_r_group_number(attachment_atom2) == no_r_group_num);
    BOOST_TEST(get_next_attachment_point_number(mol) == 3);

    model.addAtom(Element::N, RDGeom::Point3D(5.0, 6.0, 0.0));
    const auto* n_atom = mol->getAtomWithIdx(3);

    model.addAttachmentPoint(RDGeom::Point3D(7.0, 8.0, 0.0), n_atom);
    const auto* attachment_atom3 = mol->getAtomWithIdx(4);
    BOOST_TEST(attachment_atom3->getAtomicNum() == 0);
    BOOST_TEST(is_attachment_point(attachment_atom3));
    BOOST_TEST(get_attachment_point_number(attachment_atom3) == 3);
    BOOST_TEST(get_next_attachment_point_number(mol) == 4);

    n_atom = mol->getAtomWithIdx(3);
    model.addAttachmentPoint(RDGeom::Point3D(9.0, 0.0, 0.0), n_atom);
    const auto* attachment_atom4 = mol->getAtomWithIdx(5);
    BOOST_TEST(attachment_atom4->getAtomicNum() == 0);
    BOOST_TEST(is_attachment_point(attachment_atom4));
    BOOST_TEST(get_attachment_point_number(attachment_atom4) == 4);
    BOOST_TEST(get_next_attachment_point_number(mol) == 5);
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 4);

    // make sure that deleting an attachment point renumbers all remaining
    // attachment points and removes the attachment point bond
    attachment_atom2 = mol->getAtomWithIdx(2);
    model.remove({attachment_atom2}, {}, {}, {}, {});
    attachment_atom = mol->getAtomWithIdx(1);
    attachment_atom3 = mol->getAtomWithIdx(3);
    attachment_atom4 = mol->getAtomWithIdx(4);
    BOOST_TEST(get_attachment_point_number(attachment_atom) == 1);
    BOOST_TEST(get_attachment_point_number(attachment_atom3) == 2);
    BOOST_TEST(get_attachment_point_number(attachment_atom4) == 3);
    BOOST_TEST(get_next_attachment_point_number(mol) == 4);
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);

    model.remove({attachment_atom, attachment_atom3}, {}, {}, {}, {});
    attachment_atom4 = mol->getAtomWithIdx(2);
    BOOST_TEST(get_attachment_point_number(attachment_atom4) == 1);
    BOOST_TEST(get_next_attachment_point_number(mol) == 2);
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 1);

    // make sure that deleting an attachment point bond also removes the
    // attachment point itself
    auto* ap4_bond = *(mol->atomBonds(attachment_atom4).begin());
    model.remove({}, {ap4_bond}, {}, {}, {});
    BOOST_TEST(get_next_attachment_point_number(mol) == 1);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.undo();
    BOOST_TEST(get_next_attachment_point_number(mol) == 2);
    BOOST_TEST(mol->getNumAtoms() == 3);
    BOOST_TEST(mol->getNumBonds() == 1);

    undo_stack.undo();
    BOOST_TEST(get_next_attachment_point_number(mol) == 4);
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
}

BOOST_AUTO_TEST_CASE(test_removeAtom)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    c_atom = mol->getAtomWithIdx(0);
    auto* n_atom = mol->getAtomWithIdx(1);

    model.remove({c_atom}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 1);
    n_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));
    model.remove({n_atom}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    n_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(c_atom->getSymbol() == "C");
    check_coords(mol, 0, 1.0, 2.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(0)) == c_atom);
    BOOST_TEST(model.getTagForAtom(c_atom) == AtomTag(0));
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    n_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(n_atom->getSymbol() == "N");
    check_coords(mol, 0, 3.0, 4.0);
    BOOST_TEST(model.getAtomFromTag(AtomTag(1)) == n_atom);
    BOOST_TEST(model.getTagForAtom(n_atom) == AtomTag(1));
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
}

BOOST_AUTO_TEST_CASE(test_addBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    c_atom = mol->getAtomWithIdx(0);
    const auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(mol->getNumBonds() == 0);

    model.addBond(c_atom, n_atom, RDKit::Bond::BondType::SINGLE,
                  RDKit::Bond::BondDir::BEGINDASH);
    BOOST_TEST(mol->getNumBonds() == 1);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(BondTag(0)) == bond);
    BOOST_TEST(model.getTagForBond(bond) == BondTag(0));
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINDASH);

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(BondTag(0)) == bond);
    BOOST_TEST(model.getTagForBond(bond) == BondTag(0));
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINDASH);
}

BOOST_AUTO_TEST_CASE(test_addBond_query)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    const auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(mol->getNumBonds() == 0);

    auto bond_query = std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(
        RDKit::makeBondNullQuery());
    model.addBond(c_atom, n_atom, bond_query);
    BOOST_TEST(mol->getNumBonds() == 1);
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->hasQuery());
    BOOST_TEST(bond->getQuery()->getTypeLabel() == bond_query->getTypeLabel());

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->hasQuery());
    BOOST_TEST(bond->getQuery()->getTypeLabel() == bond_query->getTypeLabel());
}

BOOST_AUTO_TEST_CASE(test_addNonMolecularObject)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(!model.hasReactionArrow());
    BOOST_TEST(model.getReactionArrow() == nullptr);
    BOOST_TEST(model.getNonMolecularObjects().empty());

    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(!model.hasReactionArrow());
    BOOST_TEST(model.getReactionArrow() == nullptr);
    auto non_mol_objects = model.getNonMolecularObjects();
    BOOST_TEST(non_mol_objects.size() == 1);
    auto* plus1 = &model.m_pluses[0];
    BOOST_TEST(plus1->getType() == NonMolecularType::RXN_PLUS);
    check_coords(plus1->getCoords(), 1.0, 2.0);
    BOOST_TEST(model.m_tag_to_non_molecular_object[NonMolecularTag(0)] ==
               plus1);

    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(3.0, 4.0, 0.0));
    BOOST_TEST(!model.hasReactionArrow());
    non_mol_objects = model.getNonMolecularObjects();
    BOOST_TEST(non_mol_objects.size() == 2);
    plus1 = &model.m_pluses[0];
    BOOST_TEST(plus1->getType() == NonMolecularType::RXN_PLUS);
    check_coords(plus1->getCoords(), 1.0, 2.0);
    BOOST_TEST(model.m_tag_to_non_molecular_object[NonMolecularTag(0)] ==
               plus1);
    auto* plus2 = &model.m_pluses[1];
    BOOST_TEST(plus2->getType() == NonMolecularType::RXN_PLUS);
    check_coords(plus2->getCoords(), 3.0, 4.0);
    BOOST_TEST(model.m_tag_to_non_molecular_object[NonMolecularTag(1)] ==
               plus2);

    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(5.0, 6.0, 0.0));
    BOOST_TEST(model.hasReactionArrow());
    non_mol_objects = model.getNonMolecularObjects();
    BOOST_TEST(non_mol_objects.size() == 3);
    auto* arrow = model.getReactionArrow();
    BOOST_TEST(non_mol_objects.count(arrow) == 1);
    BOOST_TEST(arrow->getType() == NonMolecularType::RXN_ARROW);
    check_coords(arrow->getCoords(), 5.0, 6.0);
    BOOST_TEST(model.m_tag_to_non_molecular_object[NonMolecularTag(2)] ==
               arrow);

    // make sure that we can't add a second arrow
    BOOST_CHECK_THROW(
        model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                    RDGeom::Point3D(7.0, 8.0, 0.0)),
        std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_removeBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    c_atom = mol->getAtomWithIdx(0);
    const auto* n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(mol->getNumBonds() == 0);

    model.addBond(c_atom, n_atom, RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* bond = mol->getBondWithIdx(0);

    model.remove({}, {bond}, {}, {}, {});
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumBonds() == 1);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(model.getBondFromTag(BondTag(0)) == bond);
    BOOST_TEST(model.getTagForBond(bond) == BondTag(0));
    c_atom = mol->getAtomWithIdx(0);
    n_atom = mol->getAtomWithIdx(1);
    BOOST_TEST(bond->getBeginAtom() == c_atom);
    BOOST_TEST(bond->getEndAtom() == n_atom);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);

    undo_stack.redo();
    BOOST_TEST(mol->getNumBonds() == 0);
}

BOOST_AUTO_TEST_CASE(test_removeNonMolecularObject)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(5.0, 6.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(3.0, 4.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(7.0, 8.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 3);
    BOOST_TEST(model.hasReactionArrow());

    auto* arrow = model.getReactionArrow();
    model.remove({}, {}, {}, {}, {arrow});
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 2);
    BOOST_TEST(!model.hasReactionArrow());

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 3);
    BOOST_TEST(model.hasReactionArrow());

    auto* plus2 = &model.m_pluses[1];
    model.remove({}, {}, {}, {}, {plus2});
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 2);
    BOOST_TEST(model.hasReactionArrow());

    auto* plus1 = &model.m_pluses[0];
    model.remove({}, {}, {}, {}, {plus1});
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 1);
    BOOST_TEST(model.hasReactionArrow());

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 2);
    BOOST_TEST(model.hasReactionArrow());

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getNonMolecularObjects().size() == 3);
    BOOST_TEST(model.hasReactionArrow());
}

/**
 * Ensure that removing atoms and bonds in the same call works as expected
 */
BOOST_AUTO_TEST_CASE(test_removeAtomsAndBonds)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CN"));
    model.addMol(*mol_to_add);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);

    const auto* c_atom = mol->getAtomWithIdx(0);
    const auto* n_atom = mol->getAtomWithIdx(1);
    const auto* bond = mol->getBondWithIdx(0);
    model.remove({c_atom, n_atom}, {bond}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
}

/**
 * Make sure that we remove a variable attachment bond and its dummy atom
 * whenever the bond is invalidated
 */
BOOST_AUTO_TEST_CASE(test_remove_variable_attachment_bond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, VAR_ATTACH_MOL);
    const RDKit::ROMol* mol = model.getMol();

    auto sanity_check = [&mol]() {
        BOOST_TEST(is_variable_attachment_bond(mol->getBondWithIdx(6)));
        BOOST_TEST(mol->getNumAtoms() == 8);
        BOOST_TEST(mol->getNumBonds() == 7);
    };
    auto test_no_variable_attachment_bond_nor_dummy_atom = [&mol]() {
        BOOST_TEST(std::none_of(
            mol->bonds().begin(), mol->bonds().end(),
            [](auto* bond) { return is_variable_attachment_bond(bond); }));
        BOOST_TEST(
            std::none_of(mol->atoms().begin(), mol->atoms().end(),
                         [](auto* atom) { return atom->getAtomicNum() == 0; }));
    };

    sanity_check();

    // delete the variable attachment bond itself
    model.remove({}, {mol->getBondWithIdx(6)}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 7);
    BOOST_TEST(mol->getNumBonds() == 6);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete one of the variable attachment atoms
    model.remove({mol->getAtomWithIdx(2)}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 4);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete a bond between two variable attachment atoms
    model.remove({}, {mol->getBondWithIdx(1)}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 7);
    BOOST_TEST(mol->getNumBonds() == 5);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete a bond involving one variable attachment atoms
    model.remove({}, {mol->getBondWithIdx(0)}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 7);
    BOOST_TEST(mol->getNumBonds() == 5);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete an atom bound to one of the variable attachment atoms, which
    // implicitly deletes a bond from a variable attachment atom
    model.remove({mol->getAtomWithIdx(0)}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 4);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete the variable attachment bond's non-dummy atom, which implicitly
    // deletes the variable attachment bond
    model.remove({mol->getAtomWithIdx(6)}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 6);
    test_no_variable_attachment_bond_nor_dummy_atom();
    undo_stack.undo();
    sanity_check();

    // delete an atom that has nothing to do with the variable attachment bond
    // and make sure that the variable attachment bond isn't deleted
    model.remove({mol->getAtomWithIdx(5)}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 7);
    BOOST_TEST(mol->getNumBonds() == 5);
    BOOST_TEST(is_variable_attachment_bond(mol->getBondWithIdx(4)));
    undo_stack.undo();
    sanity_check();

    // delete a bond that has nothing to do with the variable attachment bond
    // and make sure that the variable attachment bond isn't deleted
    model.remove({}, {mol->getBondWithIdx(4)}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(mol->getNumBonds() == 6);
    BOOST_TEST(is_variable_attachment_bond(mol->getBondWithIdx(5)));
    undo_stack.undo();
    sanity_check();
}

BOOST_AUTO_TEST_CASE(test_addMol)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
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
        BOOST_TEST(model.getAtomFromTag(AtomTag(i)) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == AtomTag(i));
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(BondTag(i)) == bond);
        BOOST_TEST(model.getTagForBond(bond) == BondTag(i));
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
        BOOST_TEST(model.getAtomFromTag(AtomTag(i)) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == AtomTag(i));
    }
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(BondTag(i)) == bond);
        BOOST_TEST(model.getTagForBond(bond) == BondTag(i));
    }
}

BOOST_AUTO_TEST_CASE(test_addMol_attachment_points)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    import_mol_text(&model, "C*.N* |$;_AP7;;_AP3$|");
    BOOST_TEST(mol->getNumAtoms() == 4);
    BOOST_TEST(mol->getNumBonds() == 2);
    auto* ap7_atom = mol->getAtomWithIdx(1);
    auto* ap3_atom = mol->getAtomWithIdx(3);
    // addMol should automatically renumber the attachment points
    BOOST_TEST(get_attachment_point_number(ap7_atom) == 2);
    BOOST_TEST(get_attachment_point_number(ap3_atom) == 1);

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 4);
    BOOST_TEST(mol->getNumBonds() == 2);
    ap7_atom = mol->getAtomWithIdx(1);
    ap3_atom = mol->getAtomWithIdx(3);
    BOOST_TEST(get_attachment_point_number(ap7_atom) == 2);
    BOOST_TEST(get_attachment_point_number(ap3_atom) == 1);

    import_mol_text(&model, "C*.N* |$;_AP1;;_AP2$|");
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(mol->getNumBonds() == 4);
    ap7_atom = mol->getAtomWithIdx(1);
    ap3_atom = mol->getAtomWithIdx(3);
    auto* ap1_atom = mol->getAtomWithIdx(5);
    auto* ap2_atom = mol->getAtomWithIdx(7);
    // the existing _AP1 should be numbered before the new _AP1.  Same with _AP2
    BOOST_TEST(get_attachment_point_number(ap3_atom) == 1);
    BOOST_TEST(get_attachment_point_number(ap1_atom) == 2);
    BOOST_TEST(get_attachment_point_number(ap7_atom) == 3);
    BOOST_TEST(get_attachment_point_number(ap2_atom) == 4);
}

/**
 * Make sure that we can load molecules of up to MAX_NUM_ATOMS_FOR_IMPORT size,
 * but not larger
 */
BOOST_AUTO_TEST_CASE(test_addMol_size_cutoff)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, std::string(MAX_NUM_ATOMS_FOR_IMPORT, 'C'));
    BOOST_CHECK_THROW(
        import_mol_text(&model, std::string(MAX_NUM_ATOMS_FOR_IMPORT + 1, 'C')),
        std::runtime_error);
}

/**
 * Make sure that the selection is updated when a molecule is added
 */
BOOST_AUTO_TEST_CASE(test_addMol_move_selection)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCC"));
    model.addMol(*mol_to_add);
    BOOST_TEST(model.getSelectedAtoms().size() == 0);
    model.addMol(*mol_to_add, "test", true, WhatChanged::MOLECULE);
    // there was no selection, so the new atoms should not be selected
    BOOST_TEST(model.getSelectedAtoms().size() == 0);

    // now select an atom and add the molecule again
    model.select({model.getAtomFromTag(AtomTag(0))}, {}, {}, {}, {},
                 SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms().size() == 1);
    model.addMol(*mol_to_add, "test", true, WhatChanged::MOLECULE);
    // the new atoms should be selected
    BOOST_TEST(model.getSelectedAtoms().size() == 4);

    // test undo/redo
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().size() == 1);
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().size() == 4);
}

BOOST_AUTO_TEST_CASE(test_clear)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    model.clear();
    BOOST_TEST(!undo_stack.count());

    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCC"));
    model.addMol(*mol_to_add);
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(7.0, 8.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(5.0, 6.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
    BOOST_TEST(model.getNonMolecularObjects().size() == 2);
    BOOST_TEST(!model.isEmpty());

    auto* atom1 = model.getAtomFromTag(AtomTag(1));
    auto* plus = &model.m_pluses[0];
    model.select({atom1}, {}, {}, {}, {plus}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));

    model.clear();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
    BOOST_TEST(model.getNonMolecularObjects().size() == 0);
    BOOST_TEST(model.isEmpty());
    BOOST_TEST(!model.hasSelection());

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 5);
    BOOST_TEST(mol->getNumBonds() == 3);
    BOOST_TEST(model.getNonMolecularObjects().size() == 2);
    atom1 = model.getAtomFromTag(AtomTag(1));
    plus = &model.m_pluses[0];
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    const RDKit::Atom* atom;
    for (int i = 1; i < 5; ++i) {
        atom = mol->getAtomWithIdx(i);
        BOOST_TEST(atom->getSymbol() == "C");
        BOOST_TEST(model.getAtomFromTag(AtomTag(i)) == atom);
        BOOST_TEST(model.getTagForAtom(atom) == AtomTag(i));
    }
    const RDKit::Bond* bond;
    for (int i = 0; i < 3; ++i) {
        bond = mol->getBondWithIdx(i);
        BOOST_TEST(model.getBondFromTag(BondTag(i)) == bond);
        BOOST_TEST(model.getTagForBond(bond) == BondTag(i));
    }

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
    BOOST_TEST(model.getNonMolecularObjects().size() == 0);
    BOOST_TEST(model.isEmpty());
    BOOST_TEST(!model.hasSelection());

    // make sure we can add stuff after a clear SKETCH-2027
    model.addAtom(Element::N, RDGeom::Point3D(3.0, 4.0, 0.0));
}

BOOST_AUTO_TEST_CASE(test_selection)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(!model.hasSelection());
    int selection_changed_count = 0;

    QSignalSpy selection_changed_spy(&model, &MolModel::selectionChanged);
    auto test_selection_changed_emitted = [&selection_changed_count,
                                           &selection_changed_spy](
                                              bool emitted, int increment = 1) {
        if (emitted)
            selection_changed_count += increment;
        BOOST_TEST(selection_changed_spy.count() == selection_changed_count);
    };

    // calls that don't change anything shouldn't add a command to the undo
    // stack or emit selectionChanged signal
    model.clearSelection();
    test_selection_changed_emitted(false);
    BOOST_TEST(!undo_stack.count());
    model.select({}, {}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(!undo_stack.count());
    model.select({}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(!undo_stack.count());

    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCCC"));
    model.addMol(*mol_to_add);
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(7.0, 8.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(5.0, 6.0, 0.0));
    BOOST_TEST(undo_stack.count() == 3);

    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(!model.hasSelection());

    // calls that don't change anything shouldn't add a command to the undo
    // stack
    model.clearSelection();
    test_selection_changed_emitted(false);
    model.select({}, {}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(undo_stack.count() == 3);
    model.select({}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(undo_stack.count() == 3);

    auto* atom1 = model.getAtomFromTag(AtomTag(1));
    auto* atom2 = model.getAtomFromTag(AtomTag(2));
    auto* bond1 = model.getBondFromTag(BondTag(1));
    auto* arrow = &model.m_arrow.value();
    auto* plus = &model.m_pluses[0];
    model.select({atom1, atom2}, {bond1}, {}, {}, {arrow}, SelectMode::SELECT);
    test_selection_changed_emitted(true); // 1
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    BOOST_TEST(model.hasSelection());

    model.select({atom1}, {bond1}, {}, {}, {arrow}, SelectMode::DESELECT);
    test_selection_changed_emitted(true); // 2
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(model.hasSelection());
    undo_stack.undo();
    test_selection_changed_emitted(true); // 3

    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    BOOST_TEST(model.hasSelection());
    undo_stack.undo();
    test_selection_changed_emitted(true); // 4
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.redo();
    test_selection_changed_emitted(true); // 5
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    BOOST_TEST(model.hasSelection());

    model.clearSelection();
    test_selection_changed_emitted(true); // 6
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.undo();
    test_selection_changed_emitted(true); // 7
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    BOOST_TEST(model.hasSelection());
    undo_stack.redo();
    test_selection_changed_emitted(true); // 8
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    BOOST_TEST(!model.hasSelection());

    model.select({atom2}, {bond1}, {}, {}, {plus}, SelectMode::SELECT);
    test_selection_changed_emitted(true); // 9
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.undo();
    test_selection_changed_emitted(true); // 10
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    undo_stack.redo();
    test_selection_changed_emitted(true); // 11
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));

    // make sure that selecting atoms/bonds that are already selected can be
    // undone correctly
    auto* bond2 = model.getBondFromTag(BondTag(2));
    model.select({atom1, atom2}, {bond1, bond2}, {}, {}, {arrow, plus},
                 SelectMode::SELECT);
    test_selection_changed_emitted(true); // 12
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() ==
               nmo_set({arrow, plus}));
    undo_stack.undo();
    test_selection_changed_emitted(true); // 13
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.redo();
    test_selection_changed_emitted(true); // 14
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() ==
               nmo_set({arrow, plus}));

    // make sure that deselecting atoms/bonds that are already deselected can be
    // undone correctly
    model.select({}, {}, {}, {}, {arrow}, SelectMode::DESELECT);
    test_selection_changed_emitted(true); // 15
    auto* atom3 = model.getAtomFromTag(AtomTag(3));
    auto* bond3 = model.getBondFromTag(BondTag(3));
    model.select({atom2, atom3}, {bond2, bond3}, {}, {}, {arrow},
                 SelectMode::DESELECT);
    test_selection_changed_emitted(true); // 16
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.undo();
    test_selection_changed_emitted(true); // 17
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1, bond2}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.redo();
    test_selection_changed_emitted(true); // 18
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));

    // toggle the selection
    model.select({atom1}, {}, {}, {}, {arrow}, SelectMode::TOGGLE);
    // toggle emits the signal twice
    test_selection_changed_emitted(true, 2); // 20
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() ==
               nmo_set({arrow, plus}));
    model.select({atom1}, {}, {}, {}, {arrow}, SelectMode::TOGGLE);
    test_selection_changed_emitted(true, 2); // 22
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.undo();
    test_selection_changed_emitted(true, 2); // 24
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() ==
               nmo_set({arrow, plus}));
    undo_stack.undo();
    test_selection_changed_emitted(true, 2); // 26
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));

    // select-only
    model.select({atom2}, {bond3}, {}, {}, {arrow}, SelectMode::SELECT_ONLY);
    // select_only emits the signal twice
    test_selection_changed_emitted(true, 2); // 28
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    undo_stack.undo();
    test_selection_changed_emitted(true, 2); // 30
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.redo();
    test_selection_changed_emitted(true, 2); // 32
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    // select-only with no atoms or bonds specified should clear the selection
    model.select({}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    // in this case only one signal is emitted
    test_selection_changed_emitted(true); // 33
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    undo_stack.undo();
    test_selection_changed_emitted(true); // 34
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond3}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
}

BOOST_AUTO_TEST_CASE(test_select_all_and_invert)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCCCC"));
    model.addMol(*mol_to_add);
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(7.0, 8.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(5.0, 6.0, 0.0));

    model.selectAll();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    BOOST_TEST(model.getSelectedNonMolecularObjects().size() == 2);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    BOOST_TEST(model.getSelectedNonMolecularObjects().size() == 2);

    auto* atom1 = model.getAtomFromTag(AtomTag(1));
    auto* atom2 = model.getAtomFromTag(AtomTag(2));
    auto* bond1 = model.getBondFromTag(BondTag(1));
    auto* arrow = &model.m_arrow.value();
    model.select({atom1, atom2}, {bond1}, {}, {}, {arrow},
                 SelectMode::SELECT_ONLY);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    model.selectAll();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    BOOST_TEST(model.getSelectedNonMolecularObjects().size() == 2);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms().size() == 5);
    BOOST_TEST(model.getSelectedBonds().size() == 4);
    BOOST_TEST(model.getSelectedNonMolecularObjects().size() == 2);
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    auto* atom3 = model.getAtomFromTag(AtomTag(3));
    auto* atom4 = model.getAtomFromTag(AtomTag(4));
    auto* atom0 = model.getAtomFromTag(AtomTag(0));
    auto* bond2 = model.getBondFromTag(BondTag(2));
    auto* bond3 = model.getBondFromTag(BondTag(3));
    auto* bond0 = model.getBondFromTag(BondTag(0));
    auto* plus = &model.m_pluses[0];
    model.invertSelection();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom3, atom4, atom0}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond2, bond3, bond0}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
    undo_stack.undo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
    undo_stack.redo();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom3, atom4, atom0}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond2, bond3, bond0}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({plus}));
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
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW,
                                RDGeom::Point3D(7.0, 8.0, 0.0));
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS,
                                RDGeom::Point3D(5.0, 6.0, 0.0));
    auto* atom1 = model.getAtomFromTag(AtomTag(1));
    auto* atom2 = model.getAtomFromTag(AtomTag(2));
    auto* bond1 = model.getBondFromTag(BondTag(1));
    auto* arrow = &model.m_arrow.value();
    model.select({atom1, atom2}, {bond1}, {}, {}, {arrow}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    model.remove({}, {bond1}, {}, {}, {arrow});
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));

    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());

    undo_stack.undo();
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));
    bond1 = model.getBondFromTag(BondTag(1));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    // removing an atom removes all of its bonds, so this removeAtomsAndBonds
    // call will implicitly remove and deselect bond1
    model.remove({model.getAtomFromTag(AtomTag(1))}, {}, {}, {}, {});
    atom2 = model.getAtomFromTag(AtomTag(2));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom2}));
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    undo_stack.undo();
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));
    bond1 = model.getBondFromTag(BondTag(1));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    model.clear();
    BOOST_TEST(model.getSelectedAtoms().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(model.getSelectedNonMolecularObjects().empty());

    undo_stack.undo();
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));
    bond1 = model.getBondFromTag(BondTag(1));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    // removing non-selected atoms and bonds shouldn't have any effect on the
    // selection
    auto* plus = &model.m_pluses[0];
    model.remove({model.getAtomFromTag(AtomTag(3))}, {}, {}, {}, {plus});
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));
    bond1 = model.getBondFromTag(BondTag(1));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));

    undo_stack.undo();
    model.remove({}, {model.getBondFromTag(BondTag(2))}, {}, {}, {});
    atom1 = model.getAtomFromTag(AtomTag(1));
    atom2 = model.getAtomFromTag(AtomTag(2));
    bond1 = model.getBondFromTag(BondTag(1));
    arrow = &model.m_arrow.value();
    BOOST_TEST(model.getSelectedAtoms() == atom_set({atom1, atom2}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({bond1}));
    BOOST_TEST(model.getSelectedNonMolecularObjects() == nmo_set({arrow}));
}

/**
 * Make sure that attachment points and their bonds are automatically selected
 * and deselected together.
 */
BOOST_AUTO_TEST_CASE(test_select_attachment_point)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CC* |$;;_AP1$|");
    const RDKit::ROMol* mol = model.getMol();
    const auto* c_atom = mol->getAtomWithIdx(0);
    const auto* ap_atom = mol->getAtomWithIdx(2);
    const auto* ap_bond = mol->getBondWithIdx(1);
    model.select({c_atom, ap_atom}, {}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({c_atom, ap_atom}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({ap_bond}));
    model.select({}, {ap_bond}, {}, {}, {}, SelectMode::DESELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({c_atom}));
    BOOST_TEST(model.getSelectedBonds().empty());

    model.select({}, {ap_bond}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({c_atom, ap_atom}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({ap_bond}));
    model.select({ap_atom}, {}, {}, {}, {}, SelectMode::TOGGLE);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({c_atom}));
    BOOST_TEST(model.getSelectedBonds().empty());

    model.select({ap_atom}, {ap_bond}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms() == atom_set({c_atom, ap_atom}));
    BOOST_TEST(model.getSelectedBonds() == bond_set({ap_bond}));
}

BOOST_AUTO_TEST_CASE(test_mutateAtoms)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    BOOST_TEST(mol->getNumAtoms() == 1);
    const auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(!c_atom->hasQuery());
    BOOST_TEST(c_atom->getSymbol() == "C");

    model.mutateAtoms({c_atom}, Element::N);
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "N");

    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "C");

    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 1);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getSymbol() == "N");

    // mutate atom to query atom
    model.mutateAtoms({c_atom}, AtomQuery::A);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->hasQuery());
    BOOST_TEST(c_atom->getQueryType() == "A");

    // mutate query atom to a different query
    model.mutateAtoms({c_atom}, AtomQuery::MH);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->hasQuery());
    BOOST_TEST(c_atom->getQueryType() == "MH");

    // mutate query atom to an element
    model.mutateAtoms({c_atom}, Element::C);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(!c_atom->hasQuery());
    BOOST_TEST(c_atom->getSymbol() == "C");

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->hasQuery());
    BOOST_TEST(c_atom->getQueryType() == "MH");

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->hasQuery());
    BOOST_TEST(c_atom->getQueryType() == "A");

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(!c_atom->hasQuery());
    BOOST_TEST(c_atom->getSymbol() == "N");
}

BOOST_AUTO_TEST_CASE(test_mutateAtom_r_group)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));
    const auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getAtomicNum() == 6);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(get_r_group_number(c_atom) == no_r_group_num);

    model.mutateRGroups({c_atom}, 2);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(c_atom) == 2);
    BOOST_TEST(get_all_r_group_numbers(mol) == std::vector<unsigned int>({2}));
    BOOST_TEST(get_next_r_group_numbers(mol, 3) ==
               std::vector<unsigned int>({1, 3, 4}));

    undo_stack.undo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getAtomicNum() == 6);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(get_r_group_number(c_atom) == no_r_group_num);

    undo_stack.redo();
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getAtomicNum() == 0);
    BOOST_TEST(get_r_group_number(c_atom) == 2);

    model.mutateAtoms({c_atom}, Element::C);
    c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(c_atom->getAtomicNum() == 6);
    BOOST_TEST(c_atom->getSymbol() == "C");
    BOOST_TEST(get_r_group_number(c_atom) == no_r_group_num);
}

BOOST_AUTO_TEST_CASE(test_mutateAtomsToHIsotopes)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CCC");
    const RDKit::ROMol* mol = model.getMol();
    std::unordered_set<const RDKit::Atom*> atoms = {mol->getAtomWithIdx(0),
                                                    mol->getAtomWithIdx(2)};
    auto to_atom = RDKit::Atom("H");
    to_atom.setIsotope(2);
    model.mutateAtoms(atoms, to_atom);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[2H]C[2H]");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "CCC");
    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[2H]C[2H]");

    atoms = {mol->getAtomWithIdx(0), mol->getAtomWithIdx(2)};
    to_atom.setIsotope(3);
    model.mutateAtoms(atoms, to_atom);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[3H]C[3H]");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[2H]C[2H]");
}

BOOST_AUTO_TEST_CASE(test_mutateAtomsCharge)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "C");
    const RDKit::ROMol* mol = model.getMol();
    std::unordered_set<const RDKit::Atom*> atoms = {mol->getAtomWithIdx(0)};
    model.adjustChargeOnAtoms(atoms, +1);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[CH3+]");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C");
    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[CH3+]");

    atoms = {mol->getAtomWithIdx(0)};
    model.adjustChargeOnAtoms(atoms, -2);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[CH3-]");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "[CH3+]");
}

BOOST_AUTO_TEST_CASE(test_mutateAtomsRadicalElectrons)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "C");
    const RDKit::ROMol* mol = model.getMol();
    std::unordered_set<const RDKit::Atom*> atoms = {mol->getAtomWithIdx(0)};
    model.adjustRadicalElectronsOnAtoms(atoms, +2);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) == "[CH2] |^2:0|");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) == "C");
    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) == "[CH2] |^2:0|");

    atoms = {mol->getAtomWithIdx(0)};
    model.adjustRadicalElectronsOnAtoms(atoms, -1);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) == "[CH3] |^1:0|");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) == "[CH2] |^2:0|");
}

/**
 * Make sure that we can mutate atoms to and from a query
 */
BOOST_AUTO_TEST_CASE(test_mutateAtomsQuery)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "C");
    auto* mol = model.getMol();
    auto* c_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(!c_atom->hasQuery());
    auto query_atom = RDKit::QueryAtom(7);
    query_atom.setQuery(RDKit::makeAtomFormalChargeQuery(5));
    BOOST_TEST(query_atom.hasQuery());
    BOOST_TEST(query_atom.getQuery()->getDescription() == "AtomFormalCharge");
    model.mutateAtoms({c_atom}, query_atom);

    auto* query_atom_from_model = mol->getAtomWithIdx(0);
    BOOST_TEST(query_atom_from_model->hasQuery());
    auto* query = query_atom_from_model->getQuery();
    BOOST_TEST(query->getDescription() == "AtomFormalCharge");
    BOOST_TEST(dynamic_cast<const RDKit::ATOM_EQUALS_QUERY*>(query)->getVal() ==
               5);

    undo_stack.undo();
    BOOST_TEST(!mol->getAtomWithIdx(0)->hasQuery());
    undo_stack.redo();
    BOOST_TEST(mol->getAtomWithIdx(0)->hasQuery());
}

/**
 * Make sure that mutateAtoms handles enhanced stereochemistry correctly
 */
BOOST_AUTO_TEST_CASE(test_mutateAtomsEnhancedStereo)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "N[C@H](C)C(=O)O");
    auto* mol = model.getMol();
    auto* atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               rdkit_extensions::EnhancedStereo(
                   RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
    auto to_atom = RDKit::Atom("C");

    model.mutateAtoms({atom}, to_atom,
                      EnhancedStereo(RDKit::StereoGroupType::STEREO_AND, 1));
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               rdkit_extensions::EnhancedStereo(
                   RDKit::StereoGroupType::STEREO_AND, 1));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
    undo_stack.undo();
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               rdkit_extensions::EnhancedStereo(
                   RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
    undo_stack.redo();
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               rdkit_extensions::EnhancedStereo(
                   RDKit::StereoGroupType::STEREO_AND, 1));
    BOOST_TEST(mol->getStereoGroups().size() == 1);

    model.mutateAtoms({atom}, to_atom,
                      EnhancedStereo(RDKit::StereoGroupType::STEREO_OR, 3));
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(
        get_enhanced_stereo_for_atom(atom) ==
        rdkit_extensions::EnhancedStereo(RDKit::StereoGroupType::STEREO_OR, 3));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
    undo_stack.undo();
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(get_enhanced_stereo_for_atom(atom) ==
               rdkit_extensions::EnhancedStereo(
                   RDKit::StereoGroupType::STEREO_AND, 1));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
    undo_stack.redo();
    atom = mol->getAtomWithIdx(1);
    BOOST_TEST(
        get_enhanced_stereo_for_atom(atom) ==
        rdkit_extensions::EnhancedStereo(RDKit::StereoGroupType::STEREO_OR, 3));
    BOOST_TEST(mol->getStereoGroups().size() == 1);
}

BOOST_AUTO_TEST_CASE(test_mutateBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtomChain(
        Element::C,
        {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(1.0, 0.0, 0.0)},
        nullptr, RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);

    model.mutateBonds({bond}, BondTool::DOUBLE);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::DOUBLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);
    BOOST_TEST(!bond->hasQuery());

    undo_stack.undo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);

    undo_stack.redo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::DOUBLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);

    model.mutateBonds({bond}, BondTool::SINGLE_UP);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 0);
    BOOST_TEST(bond->getEndAtomIdx() == 1);

    // use mutateBonds to flip a bond
    model.mutateBonds({bond}, BondTool::SINGLE_UP);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 1);
    BOOST_TEST(bond->getEndAtomIdx() == 0);

    // undo the flip
    undo_stack.undo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 0);
    BOOST_TEST(bond->getEndAtomIdx() == 1);

    undo_stack.undo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::DOUBLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);

    undo_stack.redo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);

    // mutate regular bond to query bond
    auto bond_any_query = std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(
        RDKit::makeBondNullQuery());
    model.mutateBonds({bond}, BondTool::ANY);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->hasQuery());
    BOOST_TEST(bond->getQuery()->getFullDescription() ==
               bond_any_query->getFullDescription());

    // mutate to a different query
    auto bond_sd_query = std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(
        RDKit::makeSingleOrDoubleBondQuery());
    model.mutateBonds({bond}, BondTool::SINGLE_OR_DOUBLE);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->hasQuery());
    BOOST_TEST(bond->getQuery()->getFullDescription() ==
               bond_sd_query->getFullDescription());

    // mutate query bond to regular bond
    model.mutateBonds({bond}, BondTool::DOUBLE);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(!bond->hasQuery());
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::DOUBLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);
}

BOOST_AUTO_TEST_CASE(test_setBondTopology)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CCCC");
    const RDKit::ROMol* mol = model.getMol();
    model.setBondTopology({mol->getBondWithIdx(0), mol->getBondWithIdx(1)},
                          schrodinger::sketcher::BondTopology::IN_RING);
    BOOST_TEST(get_mol_text(&model, Format::SMARTS) ==
               "[#6]-&@[#6]-&@[#6]-[#6]");
    model.setBondTopology({mol->getBondWithIdx(0), mol->getBondWithIdx(2)},
                          schrodinger::sketcher::BondTopology::NOT_IN_RING);
    BOOST_TEST(get_mol_text(&model, Format::SMARTS) ==
               "[#6]-&!@[#6]-&@[#6]-&!@[#6]");
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMARTS) ==
               "[#6]-&@[#6]-&@[#6]-[#6]");
    undo_stack.redo();
    model.setBondTopology({mol->getBondWithIdx(0), mol->getBondWithIdx(1),
                           mol->getBondWithIdx(2)},
                          schrodinger::sketcher::BondTopology::EITHER);
    BOOST_TEST(get_mol_text(&model, Format::SMARTS) == "[#6]-[#6]-[#6]-[#6]");
}

BOOST_AUTO_TEST_CASE(test_flipBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtomChain(
        Element::C,
        {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(1.0, 0.0, 0.0)},
        nullptr, RDKit::Bond::BondType::SINGLE,
        RDKit::Bond::BondDir::BEGINWEDGE);
    const auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 0);
    BOOST_TEST(bond->getEndAtomIdx() == 1);

    model.flipBond(bond);
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 1);
    BOOST_TEST(bond->getEndAtomIdx() == 0);

    undo_stack.undo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 0);
    BOOST_TEST(bond->getEndAtomIdx() == 1);

    undo_stack.redo();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::SINGLE);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);
    BOOST_TEST(bond->getBeginAtomIdx() == 1);
    BOOST_TEST(bond->getEndAtomIdx() == 0);
}

/**
 * make sure that flipBondStereo works as expected, turning wedges to dashes and
 * vice-versa as a single undoable operation
 */
BOOST_AUTO_TEST_CASE(test_flipBondStereo)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CCCC");
    const RDKit::ROMol* mol = model.getMol();
    model.mutateBonds({mol->getBondWithIdx(0), mol->getBondWithIdx(1)},
                      BondTool::SINGLE_DOWN);
    std::vector<RDKit::Bond::BondDir> dirs{RDKit::Bond::BondDir::BEGINDASH,
                                           RDKit::Bond::BondDir::BEGINDASH,
                                           RDKit::Bond::BondDir::NONE};
    for (unsigned int i = 0; i < dirs.size(); ++i) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == dirs.at(i));
    }
    model.flipBondStereo({mol->beginBonds(), mol->endBonds()});
    std::vector<RDKit::Bond::BondDir> flipped_dirs{
        RDKit::Bond::BondDir::BEGINWEDGE, RDKit::Bond::BondDir::BEGINWEDGE,
        RDKit::Bond::BondDir::NONE};
    for (unsigned int i = 0; i < flipped_dirs.size(); ++i) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == flipped_dirs.at(i));
    }
    undo_stack.undo();
    for (unsigned int i = 0; i < dirs.size(); ++i) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == dirs.at(i));
    }
    undo_stack.redo();
    for (unsigned int i = 0; i < flipped_dirs.size(); ++i) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == flipped_dirs.at(i));
    }
}

BOOST_AUTO_TEST_CASE(test_regenerate_coords, *utf::tolerance(0.01))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto mol = model.getMol();
    model.addAtomChain(Element::C, {{-4.11, 2.88, 0.0}, {-1.02, 5.85, 0.0}});

    // These original coordinates are much longer
    // than the standard bond lengths otherwise generated by coordgen
    BOOST_REQUIRE(mol->getNumBonds() == 1);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);

    // Explicitly call regenerateCoordinates and rescale
    model.regenerateCoordinates();
    mol = model.getMol();
    BOOST_REQUIRE(mol->getNumBonds() == 1);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) == 1.5);

    // Go back to the "bad" coordinates
    undo_stack.undo();
    mol = model.getMol();
    BOOST_REQUIRE(mol->getNumBonds() == 1);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);

    // Clear the scene
    undo_stack.undo();
    mol = model.getMol();
    BOOST_REQUIRE(mol->getNumAtoms() == 0);

    // And redo
    undo_stack.redo();
    mol = model.getMol();
    BOOST_REQUIRE(mol->getNumBonds() == 1);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);

    undo_stack.redo();
    mol = model.getMol();
    BOOST_REQUIRE(mol->getNumBonds() == 1);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) == 1.5);
}

BOOST_AUTO_TEST_CASE(test_clean_up_selection, *utf::tolerance(0.01))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto mol = model.getMol();
    model.addAtomChain(Element::C, {{-4.11, 2.88, 0.0},
                                    {-1.02, 5.85, 0.0},
                                    {0.0, 0.0, 0.0},
                                    {10.0, 1.0, 1.0},
                                    {20.0, 2.0, 2.0}});

    // These original coordinates are much longer
    // than the standard bond lengths otherwise generated by coordgen
    BOOST_REQUIRE(mol->getNumBonds() == 4);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) ==
               5.938);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 2, 3) ==
               10.099);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 4) ==
               10.099);

    // select the last atom
    model.select({mol->getAtomWithIdx(4)}, {}, {}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(model.getSelectedAtoms().size() == 1);
    BOOST_TEST(model.getSelectedBonds().size() == 0);

    // call regenerateCoordinates, make sure the first three bonds still have
    // the "bad" coordinates, but the last one have been rescaled
    model.cleanUpSelection();
    mol = model.getMol();
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) ==
               5.938);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 2, 3) ==
               10.099);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 4) == 1.5);

    // Go back to the "bad" coordinates
    undo_stack.undo();
    mol = model.getMol();
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) ==
               5.938);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 2, 3) ==
               10.099);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 4) ==
               10.099);

    // Clear the selection
    model.clearSelection();

    // clean up with no selection, nothing should change
    model.cleanUpSelection();
    mol = model.getMol();

    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               4.286);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) ==
               5.938);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 2, 3) ==
               10.099);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 4) ==
               10.099);

    // select everything
    model.invertSelection();

    model.cleanUpSelection();
    mol = model.getMol();
    // this time all bonds should be rescaled
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) == 1.5);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) == 1.5);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 2, 3) == 1.5);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 4) == 1.5);
}

/**
 * Make sure that regenerateCoords doesn't crash and that it properly generates
 * plus signs
 */
BOOST_AUTO_TEST_CASE(test_regenerate_coords_reaction)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto mol = model.getMol();

    // no reactants or products, just an arrow
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW, {0, 0, 0});
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.empty());
    BOOST_TEST(mol->getNumAtoms() == 0);

    // add a single reactant
    model.addAtomChain(Element::C, {{-10.0, 0.0, 0.0}, {-9.0, 0.0, 0.0}});
    BOOST_TEST(mol->getNumAtoms() == 2);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.empty());
    BOOST_TEST(mol->getNumAtoms() == 2);

    // add another reactant
    model.addAtomChain(Element::C, {{-8.0, 0.0, 0.0}, {-7.0, 0.0, 0.0}});
    BOOST_TEST(mol->getNumAtoms() == 4);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(mol->getNumAtoms() == 4);

    // add an extraneous plus
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS, {-5, 0, 0});
    BOOST_TEST(model.m_pluses.size() == 2);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(mol->getNumAtoms() == 4);

    // add a product
    model.addAtomChain(Element::C, {{25.0, 0.0, 0.0}, {26.0, 0.0, 0.0}});
    BOOST_TEST(mol->getNumAtoms() == 6);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(mol->getNumAtoms() == 6);

    // add another product
    model.addAtomChain(Element::C, {{28.0, 0.0, 0.0}, {29.0, 0.0, 0.0}});
    BOOST_TEST(mol->getNumAtoms() == 8);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.size() == 2);
    BOOST_TEST(mol->getNumAtoms() == 8);

    // test undo and redo
    undo_stack.undo();
    BOOST_TEST(model.m_pluses.size() == 1);
    undo_stack.redo();
    BOOST_TEST(model.m_pluses.size() == 2);

    // delete the reactants
    std::unordered_set<const RDKit::Atom*> reactant_atoms = {
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
        mol->getAtomWithIdx(3)};
    model.remove(reactant_atoms, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 4);
    model.regenerateCoordinates();
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(mol->getNumAtoms() == 4);
}

BOOST_AUTO_TEST_CASE(test_flipAroundAxis, *utf::tolerance(0.0001))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    model.addAtomChain(Element::C, {RDGeom::Point3D(0.0, 0.0, 0.0),
                                    RDGeom::Point3D(1.0, 1.0, 0.0)});

    // flip around x=y axis, coordinates remain the same
    model.flipAroundSegment(RDGeom::Point3D(0.0, 0.0, 0.0),
                            RDGeom::Point3D(1.0, 1.0, 0.0),
                            {mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)});
    BOOST_TEST(mol->getConformer().getAtomPos(0).x == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(0).y == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).x == 1.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).y == 1.0);

    // flip around x=0 axis, coordinates are mirrored vertically
    model.flipAroundSegment(RDGeom::Point3D(0.0, 0.0, 0.0),
                            RDGeom::Point3D(1.0, 0.0, 0.0),
                            {mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)});

    BOOST_TEST(mol->getConformer().getAtomPos(0).x == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(0).y == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).x == 1.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).y == -1.0);

    // flip around y=0 axis, coordinates are mirrored horizontally
    model.flipAroundSegment(RDGeom::Point3D(0.0, 0.0, 0.0),
                            RDGeom::Point3D(0.0, 1.0, 0.0),
                            {mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)});

    BOOST_TEST(mol->getConformer().getAtomPos(0).x == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(0).y == 0.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).x == -1.0);
    BOOST_TEST(mol->getConformer().getAtomPos(1).y == -1.0);
}

/**
 * Pass one pair of atoms to mergeAtoms and confirm that they are correctly
 * merged.
 */
BOOST_AUTO_TEST_CASE(test_merge_one_atom_pair, *utf::tolerance(0.001))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(
        RDKit::SmilesToMol("C1CCCCC1.C1CCNC1"));
    model.addMol(*mol_to_add);
    const RDKit::ROMol* mol = model.getMol();
    auto conf = mol->getConformer();
    const auto* n_atom = mol->getAtomWithIdx(9);
    BOOST_TEST(n_atom->getSymbol() == "N");
    const auto* atom_to_replace = mol->getAtomWithIdx(0);
    BOOST_TEST(atom_to_replace->getSymbol() == "C");
    auto exp_coords = conf.getAtomPos(0);
    model.mergeAtoms({{n_atom, atom_to_replace}});

    // make sure that the resulting connectivity is correct
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCN2(CC1)CCCC2");

    // the merged atom replaced atom 9, but we also deleted atom 0, so the
    // resulting merged atom should have index 8
    n_atom = mol->getAtomWithIdx(8);
    // make sure that the merged atoms takes its element from the first atom of
    // the input pair and the coordinates from the second atom of the input pair
    BOOST_TEST(n_atom->getSymbol() == "N");
    conf = mol->getConformer();
    auto actual_coords = conf.getAtomPos(8);
    auto dist = (actual_coords - exp_coords).length();
    BOOST_TEST(dist == 0.0);

    // confirm that undo and redo work as expected
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCCCC1.C1CCNC1");
    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCN2(CC1)CCCC2");
}

/**
 * Pass two pairs of atoms to mergeAtoms and confirm that both pairs are
 * correctly merged.
 */
BOOST_AUTO_TEST_CASE(test_merge_two_atom_pairs)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(
        RDKit::SmilesToMol("C1CCCCC1.C1CCNC1"));
    model.addMol(*mol_to_add);
    const RDKit::ROMol* mol = model.getMol();
    const auto* n_atom = mol->getAtomWithIdx(9);
    BOOST_TEST(n_atom->getSymbol() == "N");
    const auto* other_atom_to_move = mol->getAtomWithIdx(8);

    const auto* atom_to_replace = mol->getAtomWithIdx(0);
    const auto* other_atom_to_replace = mol->getAtomWithIdx(1);

    model.mergeAtoms({{n_atom, atom_to_replace},
                      {other_atom_to_move, other_atom_to_replace}});

    // make sure that the resulting connectivity is correct
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCN2CCCC2C1");

    // confirm that undo and redo work as expected
    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCCCC1.C1CCNC1");
    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "C1CCN2CCCC2C1");
}

/**
 * Make sure an SDF structure with rings is imported with correct bond lengths
 * (See SKETCH-1651 and SKETCH-2025.) and that it's centered.
 */
BOOST_AUTO_TEST_CASE(test_struc_with_rings, *utf::tolerance(0.001))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, STRUC_WITH_RINGS);
    auto mol = model.getMol();
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 0, 1) ==
               BOND_LENGTH);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 1, 2) ==
               BOND_LENGTH);
    BOOST_TEST(MolTransforms::getBondLength(mol->getConformer(), 3, 5) ==
               BOND_LENGTH);

    // Check that the molecule is centered at the origin
    auto centroid = MolTransforms::computeCentroid(mol->getConformer(),
                                                   /*ignore Hs*/ false);
    BOOST_TEST(centroid.x == 0.0);
    BOOST_TEST(centroid.y == 0.0);
    BOOST_TEST(centroid.z == 0.0);
}

BOOST_AUTO_TEST_CASE(test_updateExplictiHs)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    model.addAtom(Element::C, RDGeom::Point3D(1.0, 2.0, 0.0));

    auto mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 1);

    // add explicit Hs
    model.updateExplicitHs(ExplicitHActions::ADD);
    BOOST_TEST(mol->getNumAtoms() == 5);

    // undo the adding of Hs
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 1);

    // mutate to N
    const auto* c_atom = mol->getAtomWithIdx(0);
    model.mutateAtoms({c_atom}, Element::N);

    // add explicit Hs , should add one fewer H
    model.updateExplicitHs(ExplicitHActions::ADD);
    BOOST_TEST(mol->getNumAtoms() == 4);

    // delete an H, making it implicit, and leaving the N with only 2 Hs
    model.remove({mol->getAtomWithIdx(2)}, {}, {}, {}, {});
    BOOST_TEST(mol->getNumAtoms() == 3);

    // adding Hs should add the H we deleted
    model.updateExplicitHs(ExplicitHActions::ADD);
    BOOST_TEST(mol->getNumAtoms() == 4);

    // adding explicit Hs again should not change the number of atoms
    model.updateExplicitHs(ExplicitHActions::ADD);
    BOOST_TEST(mol->getNumAtoms() == 4);

    // remove explicit Hs
    model.updateExplicitHs(ExplicitHActions::REMOVE);
    BOOST_TEST(mol->getNumAtoms() == 1);

    // undo the removal of Hs
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 4);

    // remove explicit Hs twice, should not change the number of atoms
    model.updateExplicitHs(ExplicitHActions::REMOVE);
    BOOST_TEST(mol->getNumAtoms() == 1);
    model.updateExplicitHs(ExplicitHActions::REMOVE);
    BOOST_TEST(mol->getNumAtoms() == 1);

    // remove Hs that are selected and make sure that they're deselected
    model.updateExplicitHs(ExplicitHActions::ADD);
    model.selectAll();
    BOOST_TEST(mol->getNumAtoms() == 4);
    BOOST_TEST(model.getSelectedAtoms().size() == 4);
    model.updateExplicitHs(ExplicitHActions::REMOVE);
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(model.getSelectedAtoms().size() == 1);
}

/**
 * Test toggling explicit Hs to a single atom works as expected SKETCH-2125
 */
BOOST_AUTO_TEST_CASE(test_updateExplictiHsSingleAtom)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("C1CCCCC1"));
    model.addMol(*mol_to_add);
    auto mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 6);

    // toggling Hs should add two explicit Hs
    model.updateExplicitHs(ExplicitHActions::TOGGLE, {mol->getAtomWithIdx(0)});
    BOOST_TEST(mol->getNumAtoms() == 8);

    // toggling Hs again should remove the two explicit Hs
    model.updateExplicitHs(ExplicitHActions::TOGGLE, {mol->getAtomWithIdx(0)});
    BOOST_TEST(mol->getNumAtoms() == 6);
}

/**
 * Make sure that adding a fragment works as expected.  This method mimics
 * the behavior of DrawFragmentSceneTool.
 * @param core_smiles The SMILES for the core structure
 * @param frag_smiles The SMILES for the fragment to be added
 * @param exp_result_smiles The canonical SMILES of the expected output
 * structure
 * @param core_atom_idx If >= 0, mimic the DrawFragmentSceneTool behavior
 * when the mouse cursor is over the specified atom.
 * @param core_bond_idx If >= 0, mimic the DrawFragmentSceneTool behavior
 * when the mouse cursor is over the specified bond.
 */
void check_fragment_addition(const std::string& core_smiles,
                             const std::string& frag_smiles,
                             const std::string& exp_result_smiles,
                             int core_atom_idx = -1, int core_bond_idx = -1)
{
    BOOST_TEST_MESSAGE("Running check_fragment_addition:"
                       << "\n\tcore_smiles = " << core_smiles
                       << "\n\tfrag_smiles = " << frag_smiles
                       << "\n\texp_result_smiles = " << exp_result_smiles
                       << "\n\tcore_atom_idx = " << core_atom_idx
                       << "\n\tcore_bond_idx = " << core_bond_idx << "\n");

    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, core_smiles);
    auto* mol = model.getMol();
    auto frag = rdkit_extensions::to_rdkit(frag_smiles);
    prepare_mol(*frag);
    RDKit::Conformer frag_conf;
    const RDKit::Atom* core_atom = nullptr;
    if (core_atom_idx >= 0 && core_bond_idx >= 0) {
        BOOST_FAIL("core_atom_idx and core_bond_idx cannot both be greater "
                   "than 0");
    } else if (core_atom_idx >= 0) {
        // mimic what the DrawFragmentSceneTool does when the mouse cursor
        // is over an atom
        core_atom = mol->getAtomWithIdx(core_atom_idx);
        frag_conf = align_fragment_with_atom(*frag, core_atom);
    } else if (core_bond_idx >= 0) {
        // mimic what the DrawFragmentSceneTool does when the mouse cursor
        // is over a bond
        const auto* core_bond = mol->getBondWithIdx(core_bond_idx);
        std::tie(frag_conf, core_atom) =
            align_fragment_with_bond(*frag, core_bond);
    } else {
        // mimic what the DrawFragmentSceneTool does when the mouse cursor
        // is over empty space
        frag_conf = translate_fragment_center_to(*frag, {10, 10, 0});
    }
    frag->getConformer() = frag_conf;
    model.addFragment(*frag, core_atom);
    auto result_smiles = get_mol_text(&model, Format::SMILES);
    // clang-format off
    BOOST_TEST(result_smiles == exp_result_smiles,
               "FAILURE: expected result = " <<  exp_result_smiles
                   << "\n\tactual result = " << result_smiles
                   << "\n\tcore_smiles = " << core_smiles
                   << "\n\tfrag_smiles = " << frag_smiles
                   << "\n\tcore_atom_idx = " << core_atom_idx
                   << "\n\tcore_bond_idx = " << core_bond_idx);
    // clang-format on

    // make sure that all bonds in the structure are the appropriate length
    auto& mol_conf = mol->getConformer();
    for (auto* bond : mol->bonds()) {
        auto& begin_coords = mol_conf.getAtomPos(bond->getBeginAtomIdx());
        auto& end_coords = mol_conf.getAtomPos(bond->getEndAtomIdx());
        auto bond_dist = (begin_coords - end_coords).length();
        BOOST_TEST(bond_dist == BOND_LENGTH, tt::tolerance(0.01));
    }
}

BOOST_AUTO_TEST_CASE(test_adding_fragments)
{
    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, AMIDE_FRAG_SMILES,
                            "C1CCCCC1.CC(N)=O", -1, -1);
    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, AMIDE_FRAG_SMILES,
                            "NC(=O)C1CCCCC1", 0, -1);

    check_fragment_addition(METHYLCYCLOHEXANE_MOL_SMILES, AMIDE_FRAG_SMILES,
                            "NC(=O)CC1CCCCC1", 0, -1);
    check_fragment_addition(METHYLCYCLOHEXANE_MOL_SMILES, AMIDE_FRAG_SMILES,
                            "NC(=O)C1CCCCC1", -1, 0);

    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, CYCLOHEXANE_FRAG_SMILES,
                            "C1CCCCC1.C1CCCCC1", -1, -1);
    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, CYCLOHEXANE_FRAG_SMILES,
                            "C1CCC2(CC1)CCCCC2", 0, -1);
    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, CYCLOHEXANE_FRAG_SMILES,
                            "C1CCC2CCCCC2C1", -1, 0);

    check_fragment_addition(CYCLOHEXANE_MOL_SMILES, BENZENE_FRAG_SMILES,
                            "C1=CC2=C(C=C1)CCCC2", -1, 0);
    check_fragment_addition(BENZENE_MOL_SMILES, BENZENE_FRAG_SMILES,
                            "C1=CC2=C(C=C1)C=CC=C2", -1, 0);
    check_fragment_addition(BENZENE_MOL_SMILES, BENZENE_FRAG_SMILES,
                            TWO_RINGS_SMILES, -1, 1);

    // add a third ring and make sure that we don't introduce valence errors
    check_fragment_addition(TWO_RINGS_SMILES, BENZENE_FRAG_SMILES,
                            "C1=CC2=CC=CC3=C2C(=C1)C=CC3", -1, 3);
    check_fragment_addition(TWO_RINGS_SMILES, BENZENE_FRAG_SMILES,
                            "C1=CC2=CC=CC3=C2C(=C1)CC=C3", -1, 2);

    // add a benzene ring to a sulfur atom and make sure that we keep the
    // benzene's double bond
    check_fragment_addition("CC(C)S", BENZENE_FRAG_SMILES, "CC(C)S1=CC=CC=C1",
                            3, -1);
}

/**
 * Make sure that getMolForExport returns the expected molecule and properly
 * strips all internal properties and bookmarks.
 */
BOOST_AUTO_TEST_CASE(test_getMolForExport)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CC");
    auto mol = model.getMolForExport();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(mol->getAtomBookmarks()->empty());
    BOOST_TEST(mol->getBondBookmarks()->empty());
    BOOST_TEST(!mol->getAtomWithIdx(0)->hasProp("SKETCHER_TAG"));
    BOOST_TEST(!mol->getBondWithIdx(0)->hasProp("SKETCHER_TAG"));
}

BOOST_AUTO_TEST_CASE(test_getSelectedMolForExport)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CC");

    // no selection returns an empty mol
    auto mol = model.getSelectedMolForExport();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    // partial selection
    auto atom = model.getMol()->getAtomWithIdx(0);
    model.select({atom}, {}, {}, {}, {}, SelectMode::SELECT);
    mol = model.getSelectedMolForExport();
    BOOST_TEST(mol->getNumAtoms() == 1);
    BOOST_TEST(mol->getNumBonds() == 0);
    BOOST_TEST(mol->getAtomBookmarks()->empty());
    BOOST_TEST(mol->getBondBookmarks()->empty());
    BOOST_TEST(!mol->getAtomWithIdx(0)->hasProp("SKETCHER_TAG"));

    // full selection
    model.selectAll();
    mol = model.getSelectedMolForExport();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(mol->getAtomBookmarks()->empty());
    BOOST_TEST(mol->getBondBookmarks()->empty());
    BOOST_TEST(!mol->getAtomWithIdx(0)->hasProp("SKETCHER_TAG"));
    BOOST_TEST(!mol->getBondWithIdx(0)->hasProp("SKETCHER_TAG"));
}

BOOST_AUTO_TEST_CASE(test_getReactionForExport)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto* mol = model.getMol();

    // we can't get a reaction from an empty model
    BOOST_CHECK_THROW(model.getReactionForExport(), std::runtime_error);

    // add reactants
    import_mol_text(&model, "CC");
    import_mol_text(&model, "CC");
    // we need an arrow to have a reaction
    BOOST_TEST(!model.isReactantAtom(mol->getAtomWithIdx(0)));
    BOOST_TEST(!model.isProductAtom(mol->getAtomWithIdx(0)));
    BOOST_CHECK_THROW(model.getReactionForExport(), std::runtime_error);
    model.addNonMolecularObject(NonMolecularType::RXN_ARROW, {500, 0, 0});
    BOOST_TEST(model.isReactantAtom(mol->getAtomWithIdx(0)));
    BOOST_TEST(!model.isProductAtom(mol->getAtomWithIdx(0)));
    // we need a product to have a reaction
    BOOST_CHECK_THROW(model.getReactionForExport(), std::runtime_error);

    // add product
    import_mol_text(&model, "CCCC");
    auto atom_iter = mol->atoms();
    std::vector<const RDKit::Atom*> all_atoms(atom_iter.begin(),
                                              atom_iter.end());
    std::unordered_set<const RDKit::Atom*> product_atoms(all_atoms.begin() + 4,
                                                         all_atoms.end());
    // move the product so it's to the right of the arrow
    model.translateByVector({750, 0, 0}, product_atoms);
    BOOST_TEST(!model.isReactantAtom(mol->getAtomWithIdx(4)));
    BOOST_TEST(model.isProductAtom(mol->getAtomWithIdx(4)));

    auto reaction = model.getReactionForExport();
    BOOST_TEST(reaction->getReactants().size() == 2);
    BOOST_TEST(reaction->getProducts().size() == 1);

    auto reactant = reaction->getReactants()[0];
    BOOST_TEST(reactant->getAtomBookmarks()->empty());
    BOOST_TEST(reactant->getBondBookmarks()->empty());
    BOOST_TEST(!reactant->getAtomWithIdx(0)->hasProp("SKETCHER_TAG"));
    BOOST_TEST(!reactant->getBondWithIdx(0)->hasProp("SKETCHER_TAG"));

    auto product = reaction->getProducts()[0];
    BOOST_TEST(product->getAtomBookmarks()->empty());
    BOOST_TEST(product->getBondBookmarks()->empty());
    BOOST_TEST(!product->getAtomWithIdx(0)->hasProp("SKETCHER_TAG"));
    BOOST_TEST(!product->getBondWithIdx(0)->hasProp("SKETCHER_TAG"));
}

/**
 * Make sure that addReaction allows us to add a single reaction that contains
 * up to MAX_NUM_ATOMS_FOR_IMPORT atoms
 */
BOOST_AUTO_TEST_CASE(test_addReaction)
{
    auto big_reactant = std::string(MAX_NUM_ATOMS_FOR_IMPORT / 4, 'C');
    auto big_product = std::string(MAX_NUM_ATOMS_FOR_IMPORT / 2, 'C');
    auto big_reaction = big_reactant + "." + big_reactant + ">>" + big_product;
    // sanity check in case MAX_NUM_ATOMS_FOR_IMPORT isn't divisible by 4
    BOOST_REQUIRE(big_reactant.size() * 2 + big_product.size() ==
                  MAX_NUM_ATOMS_FOR_IMPORT);

    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    // reaction contains one too many atoms
    BOOST_CHECK_THROW(import_reaction_text(&model, big_reaction + "C"),
                      std::runtime_error);
    // reaction contains the exact maximum number of atoms
    import_reaction_text(&model, big_reaction);
    // we can't have two reactions at once
    BOOST_CHECK_THROW(import_reaction_text(&model, "C.C>>CC"),
                      std::runtime_error);

    model.clear();

    const std::string RXN_BLOCK = R"($RXN

      Mrv1908  110420191440

  1  1
$MOL

  Mrv1908 11041914402D

  5  4  0  0  0  0            999 V2000
   -6.0107    1.2384    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7251    0.8259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4396    1.2384    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -6.7251    0.0009    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -5.2962    0.8259    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  RGP  3   3   1   4   2   5   3
M  END
$MOL

  Mrv1908 11041914402D

  6  5  0  0  0  0            999 V2000
   -2.0310    0.7770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7454    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4599    0.7770    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -2.7454   -0.4605    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -3.5233   -0.1755    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.3165    0.3645    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  2  5  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  RGP  3   3   1   4   2   6   3
M  END
)";
    import_reaction_text(&model, RXN_BLOCK);
}

BOOST_AUTO_TEST_CASE(test_reactions)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto* mol = model.getMol();

    // we can't get a reaction from an empty model
    BOOST_CHECK_THROW(get_reaction_text(&model, Format::SMILES),
                      std::runtime_error);

    auto smarts = "CC.CC>>CCCC";
    import_reaction_text(&model, smarts);
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(get_reaction_text(&model, Format::SMILES) == smarts);
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(model.m_pluses.size() == 0);
    BOOST_TEST(!model.m_arrow.has_value());
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(model.m_pluses.size() == 1);
    BOOST_TEST(model.m_arrow.has_value());
    BOOST_TEST(get_reaction_text(&model, Format::SMILES) == smarts);
}

/**
 * Ensure that atoms and nonmolecular objects can correctly be translated
 */
BOOST_AUTO_TEST_CASE(test_translateByVector)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CC"));
    model.addMol(*mol_to_add);
    BOOST_TEST(mol->getNumAtoms() == 2);
    auto atom1_start_position = mol->getConformer().getAtomPos(0);
    auto atom2_start_position = mol->getConformer().getAtomPos(1);
    const RDGeom::Point3D plus_position(1.0, 2.0, 0.0);
    model.addNonMolecularObject(NonMolecularType::RXN_PLUS, plus_position);

    auto plus = &model.m_pluses[0];
    BOOST_TEST(plus->getType() == NonMolecularType::RXN_PLUS);
    check_coords(plus->getCoords(), plus_position.x, plus_position.y);

    auto atom1 = mol->getAtomWithIdx(0);
    // translate the first atom and the plus sign by 10, 10
    const RDGeom::Point3D translation_vector(10.0, 10.0, 0.0);
    model.translateByVector(translation_vector, {atom1}, {plus});

    // make sure that the first atom and plus sign were translated
    auto atom1_end_position = mol->getConformer().getAtomPos(0);
    check_coords(atom1_end_position,
                 atom1_start_position.x + translation_vector.x,
                 atom1_start_position.y + translation_vector.y);
    check_coords(plus->getCoords(), plus_position.x + translation_vector.x,
                 plus_position.y + translation_vector.y);
    // make sure the second atom was not translated
    check_coords(mol->getConformer().getAtomPos(1), atom2_start_position.x,
                 atom2_start_position.y);

    // undo the translation and check that the coordinates are back to their old
    // values
    undo_stack.undo();
    check_coords(mol->getConformer().getAtomPos(0), atom1_start_position.x,
                 atom1_start_position.y);
    check_coords(plus->getCoords(), plus_position.x, plus_position.y);
}

/**
 * Given a mol, confirm that all atoms have the given chirality labels, and all
 * bonds have the given type and direction
 */
void assert_wedging_and_chiral_labels(
    const RDKit::ROMol& mol,
    const std::vector<std::string>& atom_chirality_labels,
    const std::vector<RDKit::Bond::BondType>& bond_types,
    const std::vector<RDKit::Bond::BondDir>& bond_dirs)
{
    BOOST_REQUIRE(mol.getNumAtoms() == atom_chirality_labels.size());
    for (auto atom : mol.atoms()) {
        auto actual = get_atom_chirality_label(*atom);
        BOOST_TEST(actual == atom_chirality_labels[atom->getIdx()]);
    }
    BOOST_REQUIRE(mol.getNumBonds() == bond_types.size());
    BOOST_REQUIRE(mol.getNumBonds() == bond_dirs.size());
    for (auto bond : mol.bonds()) {
        BOOST_TEST(bond->getBondType() == bond_types[bond->getIdx()]);
        BOOST_TEST(bond->getBondDir() == bond_dirs[bond->getIdx()]);
    }
};

/**
 * Tests that chiral centers lacking wedging are appropriately labeled with
 * undefined labels (ie ?), both on import of molecules, and on edit.
 */
BOOST_AUTO_TEST_CASE(test_undefined_stereo)
{
    using RDKit::Bond;
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    import_mol_text(&model, "CC(C)[C@@H](C)[C@H](C)N");
    auto mol = model.getMol();
    std::vector<std::string> atom_labels = {"", "",        "", "abs (R)",
                                            "", "abs (S)", "", ""};
    std::vector<Bond::BondType> bond_types(7, Bond::BondType::SINGLE);
    std::vector<Bond::BondDir> bond_dirs(7, Bond::BondDir::NONE);
    bond_dirs[3] = Bond::BondDir::BEGINDASH;
    bond_dirs[5] = Bond::BondDir::BEGINWEDGE;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);
    model.clear();

    // Start with the mol again, but with no wedging and 2 undefined centers
    import_mol_text(&model, "CC(C)C(C)C(C)N");
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "CC(C)C(C)C(C)N");

    mol = model.getMol();
    atom_labels = {"", "", "", "(?)", "", "(?)", "", ""};
    bond_dirs = std::vector<Bond::BondDir>(7, Bond::BondDir::NONE);
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);

    model.mutateBonds({model.getMol()->getBondWithIdx(3)},
                      BondTool::SINGLE_DOWN);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "CC(C)[C@@H](C)C(C)N");

    atom_labels[3] = "abs (R)";
    bond_dirs[3] = Bond::BondDir::BEGINDASH;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);

    model.mutateBonds({model.getMol()->getBondWithIdx(5)}, BondTool::SINGLE_UP);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) ==
               "CC(C)[C@@H](C)[C@H](C)N");

    atom_labels[5] = "abs (S)";
    bond_dirs[5] = Bond::BondDir::BEGINWEDGE;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);

    // try again in the opposite direction; start with defined centers

    // Clear the dashed bond, confirm that center goes back to undefined
    model.mutateBonds({model.getMol()->getBondWithIdx(3)}, BondTool::SINGLE);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "CC(C)C(C)[C@H](C)N");

    atom_labels[3] = "(?)";
    bond_dirs[3] = Bond::BondDir::NONE;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);

    // Clear the wedged bond, both centers should be undefined
    model.mutateBonds({model.getMol()->getBondWithIdx(5)}, BondTool::SINGLE);
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == "CC(C)C(C)C(C)N");

    atom_labels[5] = "(?)";
    bond_dirs[5] = Bond::BondDir::NONE;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);
}

BOOST_AUTO_TEST_CASE(test_attachment_point_stereo)
{
    // SKETCH-1477
    std::string molblock{R"CTAB(
  Mrv1908 02082207072D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F -0.5458 -2.1523 0 0
M  V30 2 C -0.5458 -0.6123 0 0 ATTCHPT=1
M  V30 3 N -0.5458 0.9277 0 0
M  V30 4 Cl 0.9942 -0.6123 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 4 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"};

    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, molblock);
    auto mol = model.getMol();
    BOOST_REQUIRE(mol->getNumAtoms() == 5);
    BOOST_REQUIRE(mol->getNumBonds() == 4);

    // C atom should be chiral on import and after cleanup (SKETCH-1477)
    // but not have a label since the attachment is ambiguous (SKETCH-1729).
    // Bond directions may change after cleanup (SKETCH-2147)

    std::vector<std::string> atom_labels = {"", "", "", "", ""};
    std::vector<RDKit::Bond::BondType> bond_types = {
        RDKit::Bond::BondType::SINGLE,
        RDKit::Bond::BondType::SINGLE,
        RDKit::Bond::BondType::SINGLE,
        RDKit::Bond::BondType::SINGLE,
    };
    std::vector<RDKit::Bond::BondDir> bond_dirs = {
        RDKit::Bond::BondDir::NONE,
        RDKit::Bond::BondDir::BEGINWEDGE,
        RDKit::Bond::BondDir::BEGINDASH,
        RDKit::Bond::BondDir::NONE,
    };

    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);
    model.regenerateCoordinates();

    bond_dirs[2] = RDKit::Bond::BondDir::NONE;
    assert_wedging_and_chiral_labels(*mol, atom_labels, bond_types, bond_dirs);
}

/**
 * Make sure that we can select S-groups.  Also ensure that deleting a selected
 * S-group properly deselects it, and that undoing the deletion reselects it.
 */
BOOST_AUTO_TEST_CASE(test_s_group_selection)
{
    const std::string S_GROUP_SMILES = "CCCCCCCCC |Sg:n:3,4,5:1:ht:::|";
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, S_GROUP_SMILES);
    auto mol = model.getMol();
    BOOST_TEST(model.getSelectedSGroups().empty());
    BOOST_TEST(!model.hasSelection());
    auto& s_group = RDKit::getSubstanceGroups(*mol)[0];

    // select an S-group
    model.select({}, {}, {}, {&s_group}, {}, SelectMode::SELECT);
    auto selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group) == 1);
    BOOST_TEST(model.hasSelection());

    // clear selection and undo
    model.clearSelection();
    BOOST_TEST(model.getSelectedSGroups().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.undo();
    selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group) == 1);
    BOOST_TEST(model.hasSelection());

    // remove the selected S-group and undo
    model.remove({}, {}, {}, {&s_group}, {});
    BOOST_TEST(RDKit::getSubstanceGroups(*mol).empty());
    BOOST_TEST(model.getSelectedSGroups().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.undo();
    auto& s_group2 = RDKit::getSubstanceGroups(*mol)[0];
    selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group2) == 1);
    BOOST_TEST(model.hasSelection());

    // remove an atom in the S-group (which implicitly deletes the S-group) and
    // undo
    auto* atom = mol->getAtomWithIdx(4);
    model.remove({atom}, {}, {}, {}, {});
    BOOST_TEST(RDKit::getSubstanceGroups(*mol).empty());
    BOOST_TEST(model.getSelectedSGroups().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.undo();
    auto& s_group3 = RDKit::getSubstanceGroups(*mol)[0];
    selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group3) == 1);
    BOOST_TEST(model.hasSelection());

    // remove a bond at the S-group boundary (which implicitly deletes the
    // S-group) and undo
    auto* bond = mol->getBondWithIdx(2);
    model.remove({}, {bond}, {}, {}, {});
    BOOST_TEST(RDKit::getSubstanceGroups(*mol).empty());
    BOOST_TEST(model.getSelectedSGroups().empty());
    BOOST_TEST(!model.hasSelection());
    undo_stack.undo();
    auto& s_group4 = RDKit::getSubstanceGroups(*mol)[0];
    selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group4) == 1);
    BOOST_TEST(model.hasSelection());

    // remove atoms and bonds that aren't in the S-group and make sure that the
    // S-group isn't deselected
    auto* different_atom = mol->getAtomWithIdx(8);
    auto* different_bond = mol->getBondWithIdx(7);
    model.remove({different_atom}, {different_bond}, {}, {}, {});
    auto& s_group5 = RDKit::getSubstanceGroups(*mol)[0];
    selected = model.getSelectedSGroups();
    BOOST_TEST(selected.size() == 1);
    BOOST_TEST(selected.count(&s_group5) == 1);
    BOOST_TEST(model.hasSelection());
}

/**
 * Ensure that addSGroup and modifySGroup work as expected
 */
BOOST_AUTO_TEST_CASE(test_s_group_creation_and_modification)
{
    const std::string S_GROUP_SMILES = "CCCCCCCCC |Sg:n:3,4,5:1:ht:::|";
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, S_GROUP_SMILES);
    auto mol = model.getMol();
    auto& all_s_groups = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(all_s_groups.size() == 1);
    auto& s_group = all_s_groups[0];
    BOOST_TEST(get_sgroup_type(s_group) == "SRU");
    BOOST_TEST(get_repeat_pattern_label(s_group) == "HT");
    BOOST_TEST(get_polymer_label(s_group) == "1");

    // modify an S-group that came from the SMILES input
    model.modifySGroup(&s_group, SubgroupType::SRU_POLYMER,
                       RepeatPattern::HEAD_TO_HEAD, "10,11");
    auto& all_s_groups2 = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(all_s_groups2.size() == 1);
    auto& s_group2 = all_s_groups2[0];
    BOOST_TEST(get_sgroup_type(s_group2) == "SRU");
    BOOST_TEST(get_repeat_pattern_label(s_group2) == "HH");
    BOOST_TEST(get_polymer_label(s_group) == "10,11");

    // create a new S-group
    std::unordered_set<const RDKit::Atom*> atoms = {mol->getAtomWithIdx(4)};
    model.addSGroup({mol->getAtomWithIdx(4)}, SubgroupType::COPOLYMER,
                    RepeatPattern::EITHER_UNKNOWN, "co");
    auto& all_s_groups3 = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(all_s_groups3.size() == 2);
    auto& new_s_group = all_s_groups3[1];
    BOOST_TEST(get_sgroup_type(new_s_group) == "COP");
    BOOST_TEST(get_repeat_pattern_label(new_s_group) == "EU");
    BOOST_TEST(get_polymer_label(new_s_group) == "co");

    // modify the newly created S-group
    model.modifySGroup(&new_s_group, SubgroupType::SRU_POLYMER,
                       RepeatPattern::HEAD_TO_TAIL, "999");
    auto& all_s_groups4 = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(all_s_groups4.size() == 2);
    auto& new_s_group2 = all_s_groups4[1];
    BOOST_TEST(get_sgroup_type(new_s_group2) == "SRU");
    BOOST_TEST(get_repeat_pattern_label(new_s_group2) == "HT");
    BOOST_TEST(get_polymer_label(new_s_group2) == "999");

    // undo the modifications
    undo_stack.undo();
    auto& all_s_groups5 = RDKit::getSubstanceGroups(*mol);
    BOOST_TEST(all_s_groups5.size() == 2);
    auto& new_s_group3 = all_s_groups5[1];
    BOOST_TEST(get_sgroup_type(new_s_group3) == "COP");
    BOOST_TEST(get_repeat_pattern_label(new_s_group3) == "EU");
    BOOST_TEST(get_polymer_label(new_s_group3) == "co");
}

/**
 * Make sure that translateByVector and rotateByAngle don't crash on undo/redo
 * (SKETCH-2155)
 */
BOOST_AUTO_TEST_CASE(test_translate_rotate_crash)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    std::shared_ptr<RDKit::ROMol> mol_to_add(
        RDKit::SmilesToMol("CCCCC.CCCCCCC.CCCCCCCC"));
    model.addMol(*mol_to_add);
    const RDGeom::Point3D plus_position(1.0, 2.0, 0.0);

    model.addNonMolecularObject(NonMolecularType::RXN_PLUS, plus_position);
    model.selectAll();

    unsigned int iterations = 50u;
    const RDGeom::Point3D translation_vector(2.0, 10.0, 0.0);
    auto centroid = find_centroid(*mol);
    check_coords(centroid, 0, 0);
    for (unsigned int i = 0; i < iterations; i++) {
        model.translateByVector(translation_vector, model.getSelectedAtoms(),
                                model.getSelectedNonMolecularObjects());
        centroid = find_centroid(*mol);
        check_coords(centroid, translation_vector.x * (i + 1),
                     translation_vector.y * (i + 1));

        model.rotateByAngle(1.1, centroid, model.getSelectedAtoms(),
                            model.getSelectedNonMolecularObjects());
    }
    /** make sure accessing the conformer coordinates doesn't cause a crash.
     * Also check that undo properly updates the coordinates. Note that commands
     * are not being merged together because we issued alternating translate and
     * rotate commands, so we can undo them one by one
     */
    double dummy = 0;
    for (unsigned int i = 0; i < iterations; i++) {
        // undo the rotation
        undo_stack.undo();
        centroid = find_centroid(*mol);
        check_coords(centroid, translation_vector.x * (iterations - i),
                     translation_vector.y * (iterations - i));
        // undo the translation
        undo_stack.undo();
        centroid = find_centroid(*mol);
        check_coords(centroid, translation_vector.x * (iterations - i - 1),
                     translation_vector.y * (iterations - i - 1));

        dummy += mol->getConformer().getAtomPos(2).x;
    }
    for (unsigned int i = 0; i < iterations * 2; i++) {
        undo_stack.redo();
        dummy += mol->getConformer().getAtomPos(2).x;
    }
    BOOST_TEST(dummy > 0.0);
}

BOOST_AUTO_TEST_CASE(test_addHs_hydrogen_counts)
{
    // SKETCH-2146

    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "C");

    auto atom = [&model]() { return model.getMol()->getAtomWithIdx(0); };
    auto count_hs = [&atom]() {
        bool include_node_hs = true;
        return atom()->getTotalNumHs(include_node_hs);
    };

    BOOST_REQUIRE(count_hs() == 4);

    // Set charge to alter the number of Hs
    model.setAtomCharge(atom(), 1);
    BOOST_REQUIRE(count_hs() == 3);

    // Show the Hs
    model.updateExplicitHs(ExplicitHActions::ADD, {atom()});
    BOOST_REQUIRE(count_hs() == 3);

    // Reset charge to restore Hs
    model.setAtomCharge(atom(), 0);
    BOOST_TEST(count_hs() == 4);
}

BOOST_AUTO_TEST_CASE(test_wedge_bond_replacement)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    auto smiles = "CC(C)C";
    import_mol_text(&model, smiles);

    auto bond = [&model]() { return model.getMol()->getBondWithIdx(1); };
    auto begin_atom_h_count = [&bond]() {
        return bond()->getBeginAtom()->getTotalNumHs();
    };
    auto end_atom_h_count = [&bond]() {
        return bond()->getEndAtom()->getTotalNumHs();
    };

    BOOST_REQUIRE(begin_atom_h_count() == 1);
    BOOST_REQUIRE(end_atom_h_count() == 3);

    // Add a wedge bond
    model.mutateBonds({bond()}, BondTool::SINGLE_UP);
    BOOST_REQUIRE(begin_atom_h_count() == 1);
    BOOST_REQUIRE(end_atom_h_count() == 3);

    // Make it a double bond
    model.mutateBonds({bond()}, BondTool::DOUBLE);
    BOOST_TEST(begin_atom_h_count() == 0);
    BOOST_TEST(end_atom_h_count() == 2);

    // Revert to a wedge bond
    model.mutateBonds({bond()}, BondTool::SINGLE_UP);
    BOOST_REQUIRE(begin_atom_h_count() == 1);
    BOOST_REQUIRE(end_atom_h_count() == 3);

    // Make it a triple bond now
    // note this has a valence error on the starting atom!
    model.mutateBonds({bond()}, BondTool::TRIPLE);
    BOOST_TEST(begin_atom_h_count() == 0);
    BOOST_TEST(end_atom_h_count() == 1);

    // Revert to a wedge bond again
    model.mutateBonds({bond()}, BondTool::SINGLE_UP);
    BOOST_REQUIRE(begin_atom_h_count() == 1);
    BOOST_REQUIRE(end_atom_h_count() == 3);

    // Try a zero order bond now
    model.mutateBonds({bond()}, BondTool::ZERO);
    BOOST_TEST(begin_atom_h_count() == 2);
    BOOST_TEST(end_atom_h_count() == 4);
}

/**
 * Make sure that addVariableAttachmentBond adds two atoms and a variable
 * attachment bond.  Also ensure that mutating the bond preserves the variable
 * attachment properties.
 */
BOOST_AUTO_TEST_CASE(test_addVariableAttachmentBond)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "C1CCCCC1");
    const RDKit::ROMol* mol = model.getMol();
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 6);
    std::unordered_set<const RDKit::Atom*> atoms{
        mol->getAtomWithIdx(1), mol->getAtomWithIdx(2), mol->getAtomWithIdx(3)};
    model.addVariableAttachmentBond(atoms);
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(mol->getNumBonds() == 7);
    // make sure that the new atoms and bond have tags
    BOOST_TEST(model.getTagForAtom(mol->getAtomWithIdx(6)) == 6);
    BOOST_TEST(model.getTagForAtom(mol->getAtomWithIdx(7)) == 7);
    const auto* bond = mol->getBondWithIdx(6);
    BOOST_TEST(model.getTagForBond(bond) == 6);
    BOOST_TEST(is_variable_attachment_bond(bond));
    atoms = {mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
             mol->getAtomWithIdx(3)};
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);

    // sanity check undo and redo
    undo_stack.undo();
    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 6);
    undo_stack.redo();
    BOOST_TEST(mol->getNumAtoms() == 8);
    BOOST_TEST(mol->getNumBonds() == 7);
    bond = mol->getBondWithIdx(6);
    atoms = {mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
             mol->getAtomWithIdx(3)};
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);

    // mutate the bond and make sure that it remains a variable attachment bond
    model.mutateBonds({bond}, BondTool::DOUBLE);
    bond = mol->getBondWithIdx(6);
    BOOST_TEST(is_variable_attachment_bond(bond));
    BOOST_TEST(bond->getBondType() == RDKit::Bond::BondType::DOUBLE);
    atoms = {mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
             mol->getAtomWithIdx(3)};
    BOOST_TEST(get_variable_attachment_atoms(bond) == atoms);
}
BOOST_AUTO_TEST_CASE(test_enhanced_stereo_groups_smiles)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    // Round-tripped centers through the sketcher should mark them as ABS
    import_mol_text(&model, "CC(C)[C@@H](C)[C@H](C)N");
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "CC(C)[C@@H](C)[C@H](C)N |a:3,5|");
    // make sure that the stereo groups are preserved on old and new molecules
    import_mol_text(&model, "N[C@@H](CS)C(=O)O");
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "CC(C)[C@@H](C)[C@H](C)N.N[C@@H](CS)C(=O)O |a:3,5,9|");
    model.clear();

    // Ensure that enh stereo group numbers are preserved through the sketcher
    auto smiles_with_enh_stereo =
        "C[C@H](N)[C@H](C)[C@H](C)[C@@H](N)[C@H](C)N |a:1,3,o3:7,9,&1:5|";
    import_mol_text(&model, smiles_with_enh_stereo);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               smiles_with_enh_stereo);

    // SKETCH-2191
    // Changing the last atom to carbon removes one of the atoms from OR3
    model.mutateAtoms({model.getMol()->getAtomWithIdx(11)}, Element::C);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "CC(C)[C@H](N)[C@H](C)[C@@H](C)[C@H](C)N |a:7,9,o3:3,&1:5|");
    // Changing the central wedge to double bond removes AND1 entirely
    model.mutateBonds({model.getMol()->getBondWithIdx(5)}, BondTool::DOUBLE);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "C=C([C@@H](C)[C@H](C)N)[C@H](N)C(C)C |a:2,4,o3:7|");
}

BOOST_AUTO_TEST_CASE(test_enhanced_stereo_groups_mdl)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    // Exercise import from v2000 molblocks with and without the chiral flag
    std::string v2k_mb = R"MDL(
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1
  1  3  1  0
  1  4  1  0
M  END
)MDL";
    // chiral flag off, AND center
    import_mol_text(&model, v2k_mb);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "C[C@H](O)Cl |&1:1|");
    model.clear();

    // chiral flag on, ABS center
    boost::algorithm::replace_all(v2k_mb,
                                  "  4  3  0  0  0  0  0  0  0  0999 V2000",
                                  "  4  3  0  0  1  0  0  0  0  0999 V2000");
    import_mol_text(&model, v2k_mb);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "C[C@H](O)Cl |a:1|");
    model.clear();

    // Exercise import from v3000 molblocks with and without the chiral flag
    std::string v3k_mb = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F -2.216000 -2.382857 0.000000 0
M  V30 2 C -2.218857 -0.954286 0.000000 0
M  V30 3 C -0.983143 -0.237429 0.000000 0
M  V30 4 C 0.255429 -0.949429 0.000000 0
M  V30 5 Cl 0.258286 -2.377714 0.000000 0
M  V30 6 C -3.457429 -0.242571 0.000000 0
M  V30 7 C -0.986000 1.191143 0.000000 0
M  V30 8 C 1.491143 -0.232571 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 2 6 CFG=1
M  V30 6 1 3 7 CFG=1
M  V30 7 1 4 8 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 MDLV30/STERAC1 ATOMS=(1 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)MDL";
    // chiral flag off, the unlabeled center becomes AND2
    import_mol_text(&model, v3k_mb);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "C[C@H](F)[C@H](C)[C@H](C)Cl |a:1,&2:3,&1:5|");
    model.clear();

    // chiral flag on, the unlabeled center becomes ABS
    boost::algorithm::replace_all(v3k_mb, "M  V30 COUNTS 8 7 0 0 0",
                                  "M  V30 COUNTS 8 7 0 0 1");
    import_mol_text(&model, v3k_mb);
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "C[C@H](F)[C@@H](C)[C@H](C)Cl |a:1,3,&1:5|");
}

BOOST_AUTO_TEST_CASE(test_enhanced_stereo_groups_mdl_canonicalization)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    // Showcase how SMILES canonicalization differs for enhanced stereo groups
    std::string molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N 1.299038 2.250000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 0.000000 0.000000 0.000000 0
M  V30 4 S -1.299038 0.750000 0.000000 0
M  V30 5 C 2.598076 -0.000000 0.000000 0
M  V30 6 O 2.598076 -1.500000 0.000000 0
M  V30 7 O 3.897114 0.750000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 2 5 6
M  V30 6 1 5 7
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";

    auto get_regular_smiles = [](const auto& mol_model) {
        auto mol = mol_model.getMolForExport();
        return rdkit_extensions::to_string(*mol, Format::SMILES);
    };

    // chiral flag on, ABS center
    import_mol_text(&model, molblock);
    // NOTE: standard SMILES matches extended SMILES canonicalization
    BOOST_TEST(get_regular_smiles(model) == "N[C@@H](CS)C(=O)O");
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "N[C@@H](CS)C(=O)O |a:1|");
    model.clear();

    // chiral flag off, AND center
    boost::algorithm::replace_all(molblock, "M  V30 COUNTS 7 6 0 0 1",
                                  "M  V30 COUNTS 7 6 0 0 0");
    import_mol_text(&model, molblock);
    BOOST_TEST(get_regular_smiles(model) == "N[C@@H](CS)C(=O)O");
    // NOTE: extended SMILES are canonicalized where parity is reversed!
    BOOST_TEST(get_mol_text(&model, Format::EXTENDED_SMILES) ==
               "N[C@H](CS)C(=O)O |&1:1|");
}

/**
 * test that both version of flip all works as expected, mirroring atom
 * coordinates and inverting bond dashes and wedges
 */
BOOST_AUTO_TEST_CASE(test_flip_all, *utf::tolerance(0.001))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const RDKit::ROMol* mol = model.getMol();
    std::shared_ptr<RDKit::ROMol> mol_to_add(
        RDKit::SmilesToMol("CC(C)[C@@H](C)[C@H](C)N"));
    model.addMol(*mol_to_add);
    BOOST_TEST(mol->getNumAtoms() == 8);
    auto positions = mol->getConformer().getPositions();
    std::vector<RDKit::Bond::BondDir> bond_dirs;
    for (auto* bond : mol->bonds()) {
        bond_dirs.push_back(bond->getBondDir());
    }
    std::map<RDKit::Bond::BondDir, RDKit::Bond::BondDir>
        opposite_bond_direction = {
            {RDKit::Bond::BondDir::NONE, RDKit::Bond::BondDir::NONE},
            {RDKit::Bond::BondDir::BEGINWEDGE, RDKit::Bond::BondDir::BEGINDASH},
            {RDKit::Bond::BondDir::BEGINDASH, RDKit::Bond::BondDir::BEGINWEDGE},
        };
    model.flipAllVertical();
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == atom_coords.x);
        BOOST_TEST(pos.y == -atom_coords.y);
    }
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() ==
                   opposite_bond_direction[bond_dirs[i]]);
    }
    model.flipAllVertical();
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == bond_dirs[i]);
    }
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == atom_coords.x);
        BOOST_TEST(pos.y == atom_coords.y);
    }
    model.flipAllHorizontal();
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == -atom_coords.x);
        BOOST_TEST(pos.y == atom_coords.y);
    }
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() ==
                   opposite_bond_direction[bond_dirs[i]]);
    }
}

/**
 * test that both version of flip selection work as expected, mirroring atom
 * coordinates and inverting bond dashes and wedges only for selected atoms
 */
BOOST_AUTO_TEST_CASE(test_flip_selection, *utf::tolerance(0.001))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(
        RDKit::SmilesToMol("CC(C)[C@@H](C)[C@H](C)N"));
    model.addMol(*mol_to_add);
    const RDKit::ROMol* mol = model.getMol();

    BOOST_TEST(mol->getNumAtoms() == 8);

    model.selectAll();
    // add a second molecule to make sure that only the selected atoms are
    // flipped
    model.addMol(*mol_to_add);
    mol = model.getMol();
    auto centroid = find_centroid(*mol, model.getSelectedAtoms());

    BOOST_TEST(mol->getNumAtoms() == 16);
    auto selected_atoms = model.getSelectedAtoms();
    BOOST_TEST(selected_atoms.size() == 8);
    // the selected atoms are the second half of the molecule
    for (auto* atom : selected_atoms) {
        BOOST_TEST(atom->getIdx() >= 8);
    }

    auto positions = mol->getConformer().getPositions();
    std::vector<RDKit::Bond::BondDir> bond_dirs;
    for (auto* bond : mol->bonds()) {
        bond_dirs.push_back(bond->getBondDir());
    }
    std::map<RDKit::Bond::BondDir, RDKit::Bond::BondDir>
        opposite_bond_direction = {
            {RDKit::Bond::BondDir::NONE, RDKit::Bond::BondDir::NONE},
            {RDKit::Bond::BondDir::BEGINWEDGE, RDKit::Bond::BondDir::BEGINDASH},
            {RDKit::Bond::BondDir::BEGINDASH, RDKit::Bond::BondDir::BEGINWEDGE},
        };
    model.flipSelectionHorizontal();
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        bool selected = i >= 8;
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(atom_coords.x ==
                   (selected ? 2 * centroid.x - pos.x : pos.x));
        BOOST_TEST(pos.y == atom_coords.y);
    }
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        bool selected = (mol->getBondWithIdx(i)->getBeginAtomIdx() >= 8) &&
                        (mol->getBondWithIdx(i)->getEndAtomIdx() >= 8);
        BOOST_TEST(
            mol->getBondWithIdx(i)->getBondDir() ==
            (selected ? opposite_bond_direction[bond_dirs[i]] : bond_dirs[i]));
    }
    model.flipSelectionHorizontal();
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == bond_dirs[i]);
    }
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == atom_coords.x);
        BOOST_TEST(pos.y == atom_coords.y);
    }

    model.flipSelectionVertical();
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        bool selected = i >= 8;
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == atom_coords.x);
        BOOST_TEST(atom_coords.y ==
                   (selected ? 2 * centroid.y - pos.y : pos.y));
    }
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        bool selected = (mol->getBondWithIdx(i)->getBeginAtomIdx() >= 8) &&
                        (mol->getBondWithIdx(i)->getEndAtomIdx() >= 8);
        BOOST_TEST(
            mol->getBondWithIdx(i)->getBondDir() ==
            (selected ? opposite_bond_direction[bond_dirs[i]] : bond_dirs[i]));
    }
    model.flipSelectionVertical();
    for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
        BOOST_TEST(mol->getBondWithIdx(i)->getBondDir() == bond_dirs[i]);
    }
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
        auto& pos = positions[i];
        auto atom_coords = mol->getConformer().getAtomPos(i);
        BOOST_TEST(pos.x == atom_coords.x);
        BOOST_TEST(pos.y == atom_coords.y);
    }
}

BOOST_AUTO_TEST_CASE(test_smarts_no_radicals)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    import_mol_text(&model, "[CH3]~[CH2]~*");
    BOOST_TEST(get_mol_text(&model, Format::SMARTS) == "[C&H3]~[C&H2]~*");

    auto mol = model.getMol();
    for (auto atom : mol->atoms()) {
        BOOST_TEST(atom->getNumRadicalElectrons() == 0);
    }
}

BOOST_AUTO_TEST_CASE(test_extra_bond_stereo)
{
    // CRDGEN-366

    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    const std::string molblock = R"CTAB(V2018903
  Mrv1908 02252117302D

 27 32  0  0  1  0            999 V2000
   -1.7672    2.9732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9978    3.2713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4795    2.6232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3512    2.6627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8596    2.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6425    1.2171    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3340    0.7668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3735   -0.0639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7394   -0.5722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7653   -1.3967    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0040   -1.6948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5225   -1.0468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0723   -0.3553    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6055    0.7895    0.0000 Sn  0  0  0  0  0  4  0  0  0  0  0  0
   -0.9297    1.9317    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7414    2.1489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3756    1.6404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3358    0.8097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6444    0.3595    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8615   -0.4522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3532   -1.0863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6859   -0.4780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9840    0.2914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2052    0.3198    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.9821    1.2852    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6840    2.0545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0180    0.0751    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
  9 13  1  0  0  0  0
 13 14  1  0  0  0  0
 14 15  1  0  0  0  0
  3 15  1  0  0  0  0
 15 16  1  0  0  0  0
  1 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
 19 20  1  0  0  0  0
 20 21  2  0  0  0  0
 12 21  1  0  0  0  0
 20 22  1  0  0  0  0
 22 23  2  0  0  0  0
 18 23  1  0  0  0  0
 14 24  1  0  0  0  0
  7 25  1  0  0  0  0
 25 26  2  0  0  0  0
  5 26  1  0  0  0  0
 14 27  1  0  0  0  0
M  END
$$$$

)CTAB";

    const std::unordered_map<unsigned, RDKit::Bond::BondStereo> expected = {
        {2, RDKit::Bond::BondStereo::STEREOTRANS},
        {6, RDKit::Bond::BondStereo::STEREOCIS},
        {18, RDKit::Bond::BondStereo::STEREOCIS},
        {22, RDKit::Bond::BondStereo::STEREOCIS},

    };

    import_mol_text(&model, molblock);
    BOOST_TEST(
        get_mol_text(&model, Format::SMILES) ==
        R"SMI([Cl][Sn]1([Cl])[N]2C3=CC=C2/C=C2/C=CC(=N2)/C=C2/C=C/C(=C/C4=N/C(=C\3)C=C4)[N]21)SMI");

    auto mol = model.getMol();

    unsigned num_stereo_bonds = 0;
    for (auto bond : mol->bonds()) {
        auto it = expected.find(bond->getIdx());
        if (it == expected.end()) {
            BOOST_TEST(bond->getStereo() ==
                       RDKit::Bond::BondStereo::STEREONONE);
        } else {
            BOOST_TEST(bond->getStereo() == it->second);
            ++num_stereo_bonds;
        }
    }
    BOOST_TEST(num_stereo_bonds == expected.size());
}

std::vector<std::tuple<std::string, std::string>> kekule_aromatic_smiles{
    {"C1=CC=CC=C1", "c1ccccc1"},
    {"C1=CNC=C1", "c1cc[nH]c1"},
};

BOOST_DATA_TEST_CASE(test_aromatize,
                     boost::unit_test::data::make(kekule_aromatic_smiles),
                     kekule_smiles, aromatized_smiles)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    import_mol_text(&model, kekule_smiles);
    model.aromatize();

    BOOST_TEST(get_mol_text(&model, Format::SMILES) == aromatized_smiles);

    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == kekule_smiles);

    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == aromatized_smiles);
}

BOOST_AUTO_TEST_CASE(test_kekulize)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    const std::string smiles = "c1ccccc1";
    const std::string kekulized_smiles = "C1=CC=CC=C1";

    import_mol_text(&model, smiles);
    model.kekulize();

    BOOST_TEST(get_mol_text(&model, Format::SMILES) == kekulized_smiles);

    undo_stack.undo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == smiles);

    undo_stack.redo();
    BOOST_TEST(get_mol_text(&model, Format::SMILES) == kekulized_smiles);
}

// test that when importing a molecule, a single set of valid coordinates is
// kept, or a new one is generated if none are present
BOOST_AUTO_TEST_CASE(test_mol_coords)
{
    for (unsigned int conformer_count = 0; conformer_count < 2;
         ++conformer_count) {
        QUndoStack undo_stack;
        TestMolModel model(&undo_stack);
        const RDKit::ROMol* mol = model.getMol();
        std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CC"));
        for (unsigned int i = 0; i < conformer_count; ++i) {
            RDKit::Conformer* conf = new RDKit::Conformer(2);
            conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
            conf->setAtomPos(1, RDGeom::Point3D(1, 0, 0));
            conf->set3D(false);
            mol_to_add->addConformer(conf, true);
        }
        model.addMol(*mol_to_add);
        BOOST_TEST(mol->getNumConformers() == 1);
        auto coordinates_are_zero = [](const RDKit::Conformer& conf) {
            auto positions = conf.getPositions();
            return std::all_of(
                positions.begin(), positions.end(),
                [](auto pos) { return pos.lengthSq() < 0.00001; });
        };
        BOOST_TEST(!coordinates_are_zero(mol->getConformer()));
    }
}

/**
 * Make sure that 3d conformers survive a round trip through MolModel.  (Maestro
 * sometimes uses those conformers to align the Sketcher structure to the
 * Maestro structure that it originated from.)
 */
BOOST_AUTO_TEST_CASE(test_mol_coords_3d)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCC"));
    auto conf = new RDKit::Conformer(3);
    conf->set3D(true);
    conf->setAtomPos(0, {1, 1, 1});
    conf->setAtomPos(1, {2, 2, 2});
    conf->setAtomPos(2, {3, 3, 3});
    mol_to_add->addConformer(conf, true);
    BOOST_TEST(mol_to_add->getNumConformers() == 1);
    model.addMol(*mol_to_add);
    auto new_mol = model.getMolForExport();
    BOOST_TEST(new_mol->getNumConformers() == 2);
    // make sure that the 2d conformer is first
    BOOST_TEST(!new_mol->getConformer(0).is3D());
    // make sure that the 3d conformer is still present
    BOOST_TEST(new_mol->getConformer(1).is3D());
    // make sure that the default conformer is the 2d one
    BOOST_TEST(!new_mol->getConformer().is3D());
}

/**
 * Make sure that we can generate the correct geometry after converting a single
 * bond to a triple
 */
BOOST_AUTO_TEST_CASE(test_hydridization_update, *utf::tolerance(0.05))
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CCCC");
    auto bond = model.getMol()->getBondWithIdx(1);
    model.mutateBonds({bond}, BondTool::TRIPLE);
    model.regenerateCoordinates();
    // confirm that the chain is linear
    const auto& conf = model.getMol()->getConformer();
    auto angle1 = get_angle_radians(conf.getAtomPos(0), conf.getAtomPos(1),
                                    conf.getAtomPos(2));
    BOOST_TEST(angle1 == std::numbers::pi);
    auto angle2 = get_angle_radians(conf.getAtomPos(1), conf.getAtomPos(2),
                                    conf.getAtomPos(3));
    BOOST_TEST(angle2 == std::numbers::pi);
}

/**
 * Make sure that the CIP labeler doesn't throw when it needs to examine a query
 * bond, which happens if we don't set a non-UNKNOWN bond type for the query
 * bond
 */
BOOST_AUTO_TEST_CASE(test_CIP_labeler_with_query_bonds)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "CC1CCCC[C@@H]1C");
    model.mutateBonds({model.getMol()->getBondWithIdx(0)},
                      BondTool::SINGLE_OR_DOUBLE);
    model.mutateBonds({model.getMol()->getBondWithIdx(0)}, BondTool::ANY);
    model.mutateBonds({model.getMol()->getBondWithIdx(1)},
                      BondTool::SINGLE_OR_AROMATIC);
    model.mutateBonds({model.getMol()->getBondWithIdx(1)},
                      BondTool::DOUBLE_OR_AROMATIC);
}

/**
 * Make sure that we only consider atoms to be monomeric if they came from a
 * monomeric (e.g. HELM) structure.  In particular, ensure that atomistic
 * structures read from *.mae or *.maegz are not considered monomeric, even
 * though residue information gets loaded into the mol.
 *
 * Also make sure that getMolForExport strips off the Sketcher-internal-only
 * per-atom properties and instead sets the correct mol-level HELM property.
 */
BOOST_AUTO_TEST_CASE(test_monomer_detection)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    const auto* mol = model.getMol();
    import_mol_text(&model, "CCC");
    BOOST_TEST(!contains_monomeric_atom(*mol));
    BOOST_TEST(!is_atom_monomeric(mol->getAtomWithIdx(0)));
    auto mol_for_export_1 = model.getMolForExport();
    BOOST_TEST(!rdkit_extensions::isMonomeric(*mol_for_export_1));
    BOOST_TEST(!contains_monomeric_atom(*mol_for_export_1));

    std::shared_ptr<RDKit::ROMol> mol_from_mae;
    auto mae_contents = read_testfile("1fjs_lig.mae");
    auto mae_stream = std::make_shared<std::istringstream>(mae_contents);
    auto reader = RDKit::MaeMolSupplier(mae_stream);
    mol_from_mae.reset(reader.next());
    model.addMol(*mol_from_mae);
    BOOST_TEST(!contains_monomeric_atom(*mol));
    BOOST_TEST(!is_atom_monomeric(mol->getAtomWithIdx(0)));
    BOOST_TEST(!is_atom_monomeric(mol->getAtomWithIdx(5)));
    auto mol_for_export_2 = model.getMolForExport();
    BOOST_TEST(!rdkit_extensions::isMonomeric(*mol_for_export_2));
    BOOST_TEST(!contains_monomeric_atom(*mol_for_export_2));

    import_mol_text(&model, "PEPTIDE1{A.S.D.F.G.H.W}$$$$V2.0");
    BOOST_TEST(contains_monomeric_atom(*mol));
    BOOST_TEST(!is_atom_monomeric(mol->getAtomWithIdx(0)));
    BOOST_TEST(!is_atom_monomeric(mol->getAtomWithIdx(5)));
    auto num_atoms = mol->getNumAtoms();
    BOOST_TEST(is_atom_monomeric(mol->getAtomWithIdx(num_atoms - 1)));
    auto mol_for_export_3 = model.getMolForExport();
    BOOST_TEST(rdkit_extensions::isMonomeric(*mol_for_export_3));
    // contains_monomeric_atom should return false because getMolForExport will
    // strip off the per-atom properties, since those are a Sketcher-internals
    // thing
    BOOST_TEST(!contains_monomeric_atom(*mol_for_export_3));
}

BOOST_AUTO_TEST_CASE(test_assignChiralTypesFromBondDirs_explicitHs)
{
    // Test that assignChiralTypesFromBondDirs preserves explicit H counts
    // when adding wedged/dashed bonds to cyclohexane
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    // Load cyclohexane from SMILES
    add_text_to_mol_model(model, "C1CCCCC1");
    const RDKit::ROMol* mol = model.getMol();

    BOOST_TEST(mol->getNumAtoms() == 6);
    BOOST_TEST(mol->getNumBonds() == 6);

    // Pick a carbon atom (index 0)
    auto* carbon = mol->getAtomWithIdx(0);

    // Initial state: carbon in cyclohexane has 2 ring bonds and 2 implicit
    // hydrogens
    BOOST_TEST(carbon->getDegree() == 2);
    BOOST_TEST(carbon->getNumImplicitHs() == 2);
    BOOST_TEST(carbon->getNumExplicitHs() == 0);
    BOOST_TEST(carbon->getTotalNumHs() == 2);
    BOOST_TEST(!carbon->hasValenceViolation());

    // Create a new carbon and bond it with a wedged up bond
    model.addAtom(Element::C, RDGeom::Point3D(1.5, 0.0, 0.0));
    carbon = mol->getAtomWithIdx(0);
    auto* new_carbon_1 = mol->getAtomWithIdx(6);
    model.addBond(carbon, new_carbon_1, RDKit::Bond::BondType::SINGLE,
                  RDKit::Bond::BondDir::BEGINWEDGE);

    // Assert bond was added correctly
    carbon = mol->getAtomWithIdx(0);
    auto* bond_1 = mol->getBondBetweenAtoms(0, 6);
    BOOST_REQUIRE(bond_1 != nullptr);
    BOOST_TEST(bond_1->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);

    // Verify hydrogen counts are preserved (this is what the fix addresses)
    BOOST_TEST(carbon->getDegree() == 3);
    BOOST_TEST(carbon->getNumImplicitHs() == 1);
    BOOST_TEST(carbon->getNumExplicitHs() == 0);
    BOOST_TEST(carbon->getTotalNumHs() == 1);
    BOOST_TEST(!carbon->hasValenceViolation());

    // Create a new carbon and bond it with a wedged down bond
    model.addAtom(Element::C, RDGeom::Point3D(0.0, 1.5, 0.0));
    carbon = mol->getAtomWithIdx(0);
    auto* new_carbon_2 = mol->getAtomWithIdx(7);
    model.addBond(carbon, new_carbon_2, RDKit::Bond::BondType::SINGLE,
                  RDKit::Bond::BondDir::BEGINDASH);

    // Assert bond was added correctly
    carbon = mol->getAtomWithIdx(0);
    auto* bond_2 = mol->getBondBetweenAtoms(0, 7);
    BOOST_REQUIRE(bond_2 != nullptr);
    BOOST_TEST(bond_2->getBondDir() == RDKit::Bond::BondDir::BEGINDASH);

    // Verify hydrogen counts are still preserved
    BOOST_TEST(carbon->getDegree() == 4);
    BOOST_TEST(carbon->getNumImplicitHs() == 0);
    BOOST_TEST(carbon->getNumExplicitHs() == 0);
    BOOST_TEST(carbon->getTotalNumHs() == 0);
    BOOST_TEST(!carbon->hasValenceViolation());

    // Create yet another carbon and bond it to create a pentavalent carbon
    model.addAtom(Element::C, RDGeom::Point3D(-1.5, 0.0, 0.0));
    carbon = mol->getAtomWithIdx(0);
    auto* new_carbon_3 = mol->getAtomWithIdx(8);
    model.addBond(carbon, new_carbon_3, RDKit::Bond::BondType::SINGLE,
                  RDKit::Bond::BondDir::BEGINWEDGE);

    // Assert bond was added correctly
    carbon = mol->getAtomWithIdx(0);
    auto* bond_3 = mol->getBondBetweenAtoms(0, 8);
    BOOST_REQUIRE(bond_3 != nullptr);
    BOOST_TEST(bond_3->getBondDir() == RDKit::Bond::BondDir::BEGINWEDGE);

    // Verify we now have a pentavalent carbon with valence violation
    BOOST_TEST(carbon->getDegree() == 5);
    BOOST_TEST(carbon->getNumImplicitHs() == 0);
    BOOST_TEST(carbon->getNumExplicitHs() == 0);
    BOOST_TEST(carbon->getTotalNumHs() == 0);
    BOOST_TEST(carbon->hasValenceViolation());
}

/**
 * Ensure that secondary connections can be selected and deleted independently
 * of the primary connection
 */
BOOST_AUTO_TEST_CASE(test_secondary_connections)
{
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);

    // Load cyclohexane from SMILES
    add_text_to_mol_model(model,
                          "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0");
    const RDKit::ROMol* mol = model.getMol();
    // sanity check
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    auto* bond = mol->getBondWithIdx(0);
    BOOST_TEST(contains_two_monomer_linkages(bond));

    // make sure that getAllUnselectedTags includes the tag for the secondary
    // connection
    auto unselected_tags = model.getAllUnselectedTags();
    auto unselected_bonds = std::get<1>(unselected_tags);
    BOOST_TEST(unselected_bonds.size() == 2);

    // select the primary bond
    model.select({}, {bond}, {}, {}, {}, SelectMode::SELECT);
    auto selected_bonds = model.getSelectedBonds();
    BOOST_TEST(selected_bonds.size() == 1);
    BOOST_TEST(model.getSelectedSecondaryConnections().empty());

    // select the secondary bond
    model.select({}, {}, {bond}, {}, {}, SelectMode::SELECT_ONLY);
    auto selected_secondary_connections =
        model.getSelectedSecondaryConnections();
    BOOST_TEST(selected_secondary_connections.size() == 1);
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(*selected_secondary_connections.begin() ==
               *selected_bonds.begin());

    // remove the primary bond
    model.remove({}, {bond}, {}, {}, {});
    mol = model.getMol();
    BOOST_TEST(mol->getNumBonds() == 1);
    // the selected bond is no longer a secondary connection, so it should now
    // count as a selected bond
    BOOST_TEST(model.getSelectedSecondaryConnections().empty());
    BOOST_TEST(model.getSelectedBonds().size() == 1);
    BOOST_TEST(!contains_two_monomer_linkages(mol->getBondWithIdx(0)));

    // undo the deletion
    undo_stack.undo();
    mol = model.getMol();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(contains_two_monomer_linkages(bond));

    // remove the secondary bond
    model.remove({}, {}, {bond}, {}, {});
    mol = model.getMol();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(model.getSelectedSecondaryConnections().empty());
    BOOST_TEST(model.getSelectedBonds().empty());
    BOOST_TEST(!contains_two_monomer_linkages(bond));

    // undo the deletion
    undo_stack.undo();
    mol = model.getMol();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(contains_two_monomer_linkages(bond));

    // remove both connections (thus removing the bond itself)
    model.remove({}, {bond}, {bond}, {}, {});
    mol = model.getMol();
    BOOST_TEST(mol->getNumBonds() == 0);
    BOOST_TEST(model.getSelectedSecondaryConnections().empty());
    BOOST_TEST(model.getSelectedBonds().empty());

    // undo the deletion
    undo_stack.undo();
    mol = model.getMol();
    bond = mol->getBondWithIdx(0);
    BOOST_TEST(mol->getNumBonds() == 1);
    BOOST_TEST(contains_two_monomer_linkages(bond));
}

/**
 * Test that stereo labels update in real-time when atoms are moved
 * Regression test for SKETCH-2590
 */
BOOST_AUTO_TEST_CASE(test_stereo_labels_update_on_atom_movement)
{
    // Create a molecule with a chiral center using a simpler SMILES
    // N[C@H](C)O creates a chiral center at the carbon
    QUndoStack undo_stack;
    TestMolModel model(&undo_stack);
    import_mol_text(&model, "N[C@H](C)O");
    auto* mol = model.getMol();

    // Get the chiral center atom (index 1, the carbon with N, C, and O)
    auto* chiral_atom = mol->getAtomWithIdx(1);

    // Get initial chirality label (should be "R" or "S" or "abs (R)" or "abs
    // (S)")
    std::string initial_label;
    if (chiral_atom->hasProp(RDKit::common_properties::atomNote)) {
        chiral_atom->getProp(RDKit::common_properties::atomNote, initial_label);
    }

    // The atom should have a chirality label after import
    BOOST_REQUIRE(!initial_label.empty());

    // Get the N atom (index 0) and calculate a translation to move it
    // to a position that will flip the chirality
    auto* n_atom = mol->getAtomWithIdx(0);
    auto& conf = mol->getConformer();
    RDGeom::Point3D n_pos = conf.getAtomPos(n_atom->getIdx());
    RDGeom::Point3D chiral_pos = conf.getAtomPos(chiral_atom->getIdx());

    // Move N atom to the opposite side of the chiral center (mirror position)
    // This should flip the R/S designation
    RDGeom::Point3D vector_to_move = chiral_pos - n_pos;
    vector_to_move *= 2.0; // Move to the opposite side

    // Use translateByVector to move the N atom
    std::unordered_set<const RDKit::Atom*> atoms = {n_atom};
    std::unordered_set<const NonMolecularObject*> non_mol_objs;
    model.translateByVector(vector_to_move, atoms, non_mol_objs);

    // Get the updated chirality label
    chiral_atom = mol->getAtomWithIdx(1); // Refresh pointer
    std::string updated_label;
    if (chiral_atom->hasProp(RDKit::common_properties::atomNote)) {
        chiral_atom->getProp(RDKit::common_properties::atomNote, updated_label);
    }

    // The chirality label should have been updated
    // It should either be the opposite chirality (R->S or S->R, or abs (R)->abs
    // (S)) or potentially undefined (?) if atoms are collinear
    BOOST_REQUIRE(!updated_label.empty());
    BOOST_TEST(updated_label != initial_label,
               "Chirality label should update when atom positions change, "
               "but it remained: "
                   << initial_label);
}

} // namespace sketcher
} // namespace schrodinger
