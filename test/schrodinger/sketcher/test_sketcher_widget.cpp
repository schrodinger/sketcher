#define BOOST_TEST_MODULE sketcher_widget_test

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::REACTION_TEXT_FORMATS;
using schrodinger::rdkit_extensions::TEXT_FORMATS;
using schrodinger::rdkit_extensions::to_rdkit;
using schrodinger::rdkit_extensions::to_rdkit_reaction;
using schrodinger::rdkit_extensions::to_string;

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

class TestSketcherWidget : public SketcherWidget
{
  public:
    TestSketcherWidget() : SketcherWidget(){};
    using SketcherWidget::addFromString;
    using SketcherWidget::importText;
    using SketcherWidget::m_mol_model;
    using SketcherWidget::m_scene;
    using SketcherWidget::m_sketcher_model;
    using SketcherWidget::m_undo_stack;
    using SketcherWidget::m_watermark_item;
};

BOOST_AUTO_TEST_CASE(test_addRDKitMolecule_getRDKitMolecule)
{
    TestSketcherWidget sk;
    auto mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);

    sk.addRDKitMolecule(*to_rdkit("CC"));
    mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
    sk.m_undo_stack->undo();
    mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 0);
    BOOST_TEST(mol->getNumBonds() == 0);
    sk.m_undo_stack->redo();
    mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 2);
    BOOST_TEST(mol->getNumBonds() == 1);
}

BOOST_AUTO_TEST_CASE(test_addRDKitReaction_getRDKitReaction)
{
    TestSketcherWidget sk;
    // we can't get a reaction from an empty model
    BOOST_CHECK_THROW(sk.getRDKitReaction(), std::runtime_error);

    sk.addRDKitReaction(*to_rdkit_reaction("CC.CC>>CCCC"));
    auto rxn = sk.getRDKitReaction();
    auto reactants = rxn->getReactants();
    BOOST_REQUIRE(reactants.size() == 2);
    BOOST_TEST(reactants[0]->getNumAtoms() == 2);
    BOOST_TEST(reactants[1]->getNumAtoms() == 2);
    auto products = rxn->getProducts();
    BOOST_REQUIRE(products.size() == 1);
    BOOST_TEST(products[0]->getNumAtoms() == 4);
    BOOST_TEST(sk.m_mol_model->hasReactionArrow());
    sk.m_undo_stack->undo();
    BOOST_CHECK_THROW(sk.getRDKitReaction(), std::runtime_error);
    BOOST_TEST(!sk.m_mol_model->hasReactionArrow());
    sk.m_undo_stack->redo();
    rxn = sk.getRDKitReaction();
    reactants = rxn->getReactants();
    BOOST_REQUIRE(reactants.size() == 2);
    BOOST_TEST(reactants[0]->getNumAtoms() == 2);
    BOOST_TEST(reactants[1]->getNumAtoms() == 2);
    products = rxn->getProducts();
    BOOST_REQUIRE(products.size() == 1);
    BOOST_TEST(products[0]->getNumAtoms() == 4);
    BOOST_TEST(sk.m_mol_model->hasReactionArrow());
}

BOOST_DATA_TEST_CASE(test_addFromString_getString_mol,
                     boost::unit_test::data::make(TEXT_FORMATS), format)
{
    auto mol = to_rdkit("C1=CC=CC=C1");
    auto text = to_string(*mol, format);

    TestSketcherWidget sk;
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == "C1=CC=CC=C1");
}

BOOST_DATA_TEST_CASE(test_addFromString_getString_reaction,
                     boost::unit_test::data::make(REACTION_TEXT_FORMATS),
                     format)
{
    if (format == Format::MDL_MOLV2000 || format == Format::MDL_MOLV3000) {
        // FIXME: see SKETCH-2071, returns C~C(~O)=O.C~C~O>>C~C~O~C(~C)=O
        return;
    }

    auto rxn = to_rdkit_reaction("CC(=O)O.OCC>>CC(=O)OCC");
    auto text = to_string(*rxn, format);

    TestSketcherWidget sk;
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == "CC(O)=O.CCO>>CCOC(C)=O");
}

BOOST_AUTO_TEST_CASE(test_importText_slot)
{
    TestSketcherWidget sk;
    sk.importText("c1nccc2n1ccc2", Format::SMILES);
    auto mol = sk.m_mol_model->getMol();
    BOOST_TEST_REQUIRE(mol != nullptr);
    BOOST_TEST(mol->getNumAtoms() == 9);
    unsigned num_atoms = 0;
    unsigned num_bonds = 0;
    for (auto item : sk.m_scene->items()) {
        if (item->type() == AtomItem::Type) {
            ++num_atoms;
            // make sure that this really is an AtomItem
            auto cast_item = dynamic_cast<AtomItem*>(item);
            BOOST_TEST(cast_item);
        } else if (item->type() == BondItem::Type) {
            ++num_bonds;
            // make sure that this really is a BondItem
            auto cast_item = dynamic_cast<BondItem*>(item);
            BOOST_TEST(cast_item);
        }
    }
    BOOST_TEST(num_atoms == 9);
    BOOST_TEST(num_bonds == 10);

    sk.m_mol_model->clear();
    mol = sk.m_mol_model->getMol();
    BOOST_TEST(mol->getNumAtoms() == 0);

    // import failed, exception caught, still an empty scene
    sk.importText("nonsense", Format::AUTO_DETECT);
    mol = sk.m_mol_model->getMol();
    BOOST_TEST(mol->getNumAtoms() == 0);
}

/**
 * Verify that watermark visibility is trigged on inclusion/removal of atoms
 */
BOOST_AUTO_TEST_CASE(test_watermark)
{
    TestSketcherWidget sk;
    // Without the event loop, we need to manually trigger Scene::changed
    QList<QRectF> region;
    sk.m_scene->changed(region);
    // sketcher starts out empty
    BOOST_TEST(sk.m_watermark_item->isVisible());
    sk.addFromString("c1ccncc1", Format::SMILES);
    sk.m_scene->changed(region);
    BOOST_TEST(!sk.m_watermark_item->isVisible());
    sk.m_mol_model->clear();
    sk.m_scene->changed(region);
    BOOST_TEST(sk.m_watermark_item->isVisible());
}

BOOST_AUTO_TEST_CASE(test_undo)
{
    TestSketcherWidget sk;

    // Add a structure to the scene
    sk.addFromString("[H][C@@](C)(Cl)[C@@]([H])(C)C([H])(C)Br", Format::SMILES);
    auto num_atoms = sk.m_mol_model->getMol()->getNumAtoms();
    BOOST_TEST(num_atoms > 0);

    // By default, the new molecule is not selected
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);

    // Test out `selectAll()`
    auto command_idx = sk.m_undo_stack->index();
    auto command_count = sk.m_undo_stack->count();
    sk.m_mol_model->selectAll();

    // Verify that all atoms were selected
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == num_atoms);

    // Verify that the selection command was added to the stack
    BOOST_TEST(sk.m_undo_stack->count() == command_count + 1);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx + 1);

    // Attempt to undo the selection
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.m_undo_stack->count() == command_count + 1);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);

    // Attempt to redo the selection
    sk.m_undo_stack->redo();
    BOOST_TEST(sk.m_undo_stack->count() == command_count + 1);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx + 1);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == num_atoms);

    // Manually clear the MolModel.
    sk.m_mol_model->clear();

    // Verify that the clear worked
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);

    // Verify that the clear operation was added to the stack
    BOOST_TEST(sk.m_undo_stack->count() == command_count + 2);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx + 2);

    // Attempt to undo the clear
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.m_undo_stack->count() == command_idx + 2);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx + 1);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == num_atoms);

    sk.m_undo_stack->redo();
    BOOST_TEST(sk.m_undo_stack->count() == command_count + 2);
    BOOST_TEST(sk.m_undo_stack->index() == command_idx + 2);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
}

/**
 * Verify that the Draw tool and the Erase tool get switched to the Select tool
 * if the undo stack creates a selection
 */
BOOST_AUTO_TEST_CASE(test_toolChangeOnSelection)
{
    TestSketcherWidget sk;

    // import a molecule and select all of the atoms
    sk.addFromString("CCCC", Format::SMILES);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::ATOM);
    sk.m_mol_model->selectAll();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 4);
    // Auto switches to the select tool
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);

    // undo the selection, switch to a draw tool, and then redo the selection,
    // which should automatically switch back to select
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);
    sk.m_sketcher_model->setValue(ModelKey::DRAW_TOOL, DrawTool::ATOM);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::ATOM);
    sk.m_undo_stack->redo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 4);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);

    // undo the selection, switch to the erase tool, and then redo the
    // selection, which should automatically switch back to select
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);
    sk.m_sketcher_model->setValue(ModelKey::DRAW_TOOL, DrawTool::ERASE);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::ERASE);
    sk.m_undo_stack->redo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 4);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);

    // undo the selection, switch to the move tool, and then redo the selection,
    // which should leave the move tool active
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::SELECT);
    sk.m_sketcher_model->setValue(ModelKey::DRAW_TOOL, DrawTool::MOVE_ROTATE);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::MOVE_ROTATE);
    sk.m_undo_stack->redo();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 4);
    BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::MOVE_ROTATE);
}
