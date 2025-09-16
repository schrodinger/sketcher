#define BOOST_TEST_MODULE sketcher_widget_test

#include <string>

#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::MOL_FORMATS;
using schrodinger::rdkit_extensions::RXN_FORMATS;
using schrodinger::rdkit_extensions::to_rdkit;
using schrodinger::rdkit_extensions::to_rdkit_reaction;
using schrodinger::rdkit_extensions::to_string;

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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
                     boost::unit_test::data::make(MOL_FORMATS), format)
{
    auto mol = to_rdkit("C1=CC=CC=C1");
    auto text = to_string(*mol, format);

    TestSketcherWidget sk;
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == "C1=CC=CC=C1");

    // test roundtripping stereochemistry
    std::string stereo_smiles = "CC(C)[C@@H](C)[C@H](C)N";
    sk.m_mol_model->clear();
    sk.addFromString(stereo_smiles);
    BOOST_TEST(sk.getString(Format::SMILES) == stereo_smiles);
}

BOOST_DATA_TEST_CASE(test_addFromString_getString_reaction,
                     boost::unit_test::data::make(RXN_FORMATS), format)
{
    auto rxn = to_rdkit_reaction("CC(=O)O.OCC>>CC(=O)OCC");
    auto text = to_string(*rxn, format);

    TestSketcherWidget sk;
    sk.addFromString(text);

    BOOST_TEST(sk.getString(Format::SMILES) == "CC(=O)O.CCO>>CCOC(C)=O");
}

/**
 * Make sure that molfiles are exported using 2D conformations rather than 3D
 */
BOOST_AUTO_TEST_CASE(test_2D_molfile)
{
    TestSketcherWidget sk;
    sk.addFromString("C");
    auto molfile = sk.getString(Format::MDL_MOLV3000);
    BOOST_TEST(molfile.find("2D") != std::string::npos);
    BOOST_TEST(molfile.find("3D") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(test_cut_copy_paste)
{
    TestSketcherWidget sk;
    std::string smiles{"C1=CC=CC=C1"};
    sk.importText(smiles, Format::SMILES);

    // no selection
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(!sk.m_mol_model->hasSelection());
    // nothing cut, nothing pasted
    sk.cut(Format::SMILES);
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    // nothing copied, nothing pasted
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    // everything copied, contents doubled
    sk.copy(Format::SMILES, SceneSubset::ALL);
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles + "." + smiles);
    BOOST_TEST(!sk.m_mol_model->hasSelection());
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(!sk.m_mol_model->hasSelection());

    // select a single bond; the cut/copy action should automatically select
    // atoms attached to any selected bonds
    auto bond = sk.m_mol_model->getMol()->getBondWithIdx(0);
    sk.m_mol_model->select({}, {bond}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);
    sk.cut(Format::SMILES);
    BOOST_TEST(sk.getString(Format::SMILES) == "C=CC=C");
    BOOST_TEST(!sk.m_mol_model->hasSelection());
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == "C=C.C=CC=C");
    BOOST_TEST(!sk.m_mol_model->hasSelection());
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.getString(Format::SMILES) == "C=CC=C");
    BOOST_TEST(!sk.m_mol_model->hasSelection());
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    // from the cut selection expansion
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 2);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);
    sk.m_mol_model->clearSelection();

    // repeat, but with copy; the selection should be expanded to include
    // atoms attached to any selected bonds, but the selection should persist
    bond = sk.m_mol_model->getMol()->getBondWithIdx(0);
    sk.m_mol_model->select({}, {bond}, {}, {}, SelectMode::SELECT);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 2);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles + ".C=C");
    BOOST_TEST(sk.m_mol_model->hasSelection());
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 2);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);

    // copying ALL is agnostic of an selection
    sk.copy(Format::SMILES, SceneSubset::ALL);
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(sk.m_mol_model->hasSelection());
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles + "." + smiles);
    BOOST_TEST(sk.m_mol_model->hasSelection());
    sk.m_undo_stack->undo();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 2);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 1);

    // selecting all and copying should be equivalent to doubling contents
    sk.m_mol_model->selectAll();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 6);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 6);
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.m_mol_model->hasSelection());
    sk.paste();
    BOOST_TEST(sk.getString(Format::SMILES) == smiles + "." + smiles);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 6);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 6);
    // selecting all and cutting should be equivalent to clearing the sketcher
    sk.m_mol_model->selectAll();
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 12);
    BOOST_TEST(sk.m_mol_model->getSelectedBonds().size() == 12);
    sk.cut(Format::SMILES);
    BOOST_TEST(sk.getString(Format::SMILES) == "");
    BOOST_TEST(!sk.m_mol_model->hasSelection());
}

/**
 * Make sure that pasting a string containing Windows newline characters work
 * correctly.  (RDKit parsing can't handle \r's on WASM builds, so
 * SketcherWidget must remove those manually.)
 */
BOOST_AUTO_TEST_CASE(test_paste_with_windows_newline)
{
    TestSketcherWidget sk;
    sk.setClipboardContents("CCC\n\r");
    sk.paste();
    auto mol = sk.m_mol_model->getMol();
    BOOST_TEST(mol->getNumAtoms() == 3);
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

/**
 * Verify hitting the atom chain tool doesn't change the bonds in a selection
 */
BOOST_AUTO_TEST_CASE(test_toolAtomChainTool)
{
    TestSketcherWidget sk;

    // import a molecule
    sk.addFromString("CCC", Format::SMILES);
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 0);
    auto model = sk.m_mol_model;

    // trigger the atom chain tool and make sure there's no crash (SKETCH-2113)
    sk.m_sketcher_model->pingValue(ModelKey::BOND_TOOL,
                                   QVariant::fromValue(BondTool::ATOM_CHAIN));
    // select everything
    model->selectAll();
    BOOST_TEST(model->getSelectedBonds().size() == 2);
    // trigger the atom chain tool again and make sure all bonds stay as single
    // bonds
    sk.m_sketcher_model->pingValue(ModelKey::BOND_TOOL,
                                   QVariant::fromValue(BondTool::ATOM_CHAIN));
    for (auto bond : model->getSelectedBonds()) {
        BOOST_TEST(bond->getBondType() == RDKit::Bond::SINGLE);
    }
    // trigger the double bond tool and make sure all bonds are turned into
    // double bonds
    sk.m_sketcher_model->pingValue(ModelKey::BOND_TOOL,
                                   QVariant::fromValue(BondTool::DOUBLE));

    for (auto bond : model->getSelectedBonds()) {
        BOOST_TEST(bond->getBondType() == RDKit::Bond::DOUBLE);
    }
    // trigger the atom chain tool again and make sure all bonds stay double
    sk.m_sketcher_model->pingValue(ModelKey::BOND_TOOL,
                                   QVariant::fromValue(BondTool::ATOM_CHAIN));
    for (auto bond : model->getSelectedBonds()) {
        BOOST_TEST(bond->getBondType() == RDKit::Bond::DOUBLE);
    }
}

/**
 * Verify that the select, move-rotate, and erase tools are switched to draw C
 * when the scene is emptied
 */
BOOST_AUTO_TEST_CASE(test_switch_to_C_on_empty_scene)
{
    TestSketcherWidget sk;
    for (auto drawTool :
         {DrawTool::SELECT, DrawTool::MOVE_ROTATE, DrawTool::ERASE}) {
        sk.addFromString("CCC", Format::SMILES);
        sk.m_sketcher_model->setValue(ModelKey::DRAW_TOOL, drawTool);
        BOOST_TEST(sk.m_sketcher_model->getDrawTool() == drawTool);
        sk.m_mol_model->clear();
        // Without the event loop, we need to manually trigger Scene::changed
        QList<QRectF> region;
        sk.m_scene->changed(region);
        BOOST_TEST(sk.m_sketcher_model->getDrawTool() == DrawTool::ATOM);
        BOOST_TEST(sk.m_sketcher_model->getAtomTool() == AtomTool::ELEMENT);
        BOOST_TEST(sk.m_sketcher_model->getElement() == Element::C);
    }
}

BOOST_AUTO_TEST_CASE(test_zoom_on_small_molecule)
{
    // test that zooming in on a small molecule doesn't zoom too much in
    TestSketcherWidget sk;
    sk.addFromString("C", Format::SMILES);
    auto view = sk.m_ui->view;
    view->fitToScreen();
    auto matrix = view->transform();
    float m11 = matrix.m11();
    float m22 = matrix.m22();
    BOOST_TEST(m11 <= 1.0);
    BOOST_TEST(m22 <= 1.0);
}

BOOST_DATA_TEST_CASE(test_auto_detect_through_sketcher_interface,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    std::string reference = "c1ccccc1";
    auto mol = to_rdkit(reference, Format::SMILES);
    auto text = to_string(*mol, sample);

    if (sample == Format::MAESTRO || sample == Format::INCHI ||
        sample == Format::PDB || sample == Format::XYZ) {
        // these exports force kekulization
        reference = "C1=CC=CC=C1";
    } else if (sample == Format::SMARTS || sample == Format::EXTENDED_SMARTS ||
               sample == Format::MDL_MOLV2000 ||
               sample == Format::MDL_MOLV3000 || sample == Format::MRV) {
        reference = "C1:C:C:C:C:C:1";
    }
    // all other formats should roundtrip the aromatic input

    // Check roundtripping and auto-detect
    TestSketcherWidget sk;
    sk.addFromString(text, sample);
    BOOST_TEST(sk.getString(Format::SMILES) == reference);
    sk.clear();
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == reference);
}

BOOST_DATA_TEST_CASE(test_reactions_roundtrip,
                     boost::unit_test::data::make(RXN_FORMATS))
{
    std::string smiles = "CC(=O)O.OCC>>CC(=O)OCC";
    auto reaction = to_rdkit_reaction(smiles, Format::SMILES);
    auto text = to_string(*reaction, sample);

    // Check auto-detect import through the sketcher interface as well
    TestSketcherWidget sk;
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == "CC(=O)O.CCO>>CCOC(C)=O");
}

BOOST_AUTO_TEST_CASE(testRGRoup0)
{

    TestSketcherWidget sk;

    auto molblock = read_testfile("shared_8381.sdf");

    BOOST_REQUIRE(molblock.find("RGROUPS=(1 0)") != std::string::npos);
    BOOST_REQUIRE(molblock.find("RGROUPS=(1 1)") == std::string::npos);
    BOOST_REQUIRE(molblock.find("RGROUPS=(1 6)") == std::string::npos);

    sk.addFromString(molblock, Format::MDL_MOLV3000);
    auto roundtrip_molblock = sk.getString(Format::MDL_MOLV3000);

    size_t pos{std::string::npos};
    for (int i = 2; i <= 6; ++i) {
        pos = roundtrip_molblock.find(fmt::format("RGROUPS=(1 {})", i));
        BOOST_TEST(pos != std::string::npos);
    }

    // There should be 2 R6 groups
    ++pos;
    BOOST_TEST(roundtrip_molblock.find("RGROUPS=(1 6)", pos) !=
               std::string::npos);

    BOOST_TEST(roundtrip_molblock.find("RGROUPS=(1 0)") == std::string::npos);
    BOOST_TEST(roundtrip_molblock.find("RGROUPS=(1 1)") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(testSMARTSRGRoup)
{

    TestSketcherWidget sk;

    std::string cxsmarts{
        "[H][#7](-[$([#1,*])])-[#6]-1=[#7]-[#6]-2=[#6](-[#7]-1)-"
        "[#6](Cl)=[#6]-[#6]=[#7]-2 |$;;_R1;;;;;;;;;;;$|"};

    sk.addFromString(cxsmarts, Format::SMARTS);
    BOOST_TEST(sk.getString(Format::EXTENDED_SMARTS) ==
               "[#1][#7](-[$(*)])-[#6]1=[#7]-[#6]2=[#6](-[#7]-1)-[#6](Cl)=[#6]-"
               "[#6]=[#7]-2 |$;;_R1;;;;;;;;;;$|");
}

BOOST_AUTO_TEST_CASE(testCXSMILESRGRoup)
{
    // SKETCH-1399

    TestSketcherWidget sk;

    std::string cxsmiles{"[*]C1=CC=CC=C1 |$_R1;;;;;;$,c:3,5,t:1|"};

    sk.addFromString(cxsmiles, Format::EXTENDED_SMILES);
    BOOST_TEST(sk.getString(Format::EXTENDED_SMARTS) ==
               "[#0]-[#6]1=[#6]-[#6]=[#6]-[#6]=[#6]-1 |$_R1;;;;;;$|");
}

BOOST_AUTO_TEST_CASE(testLeakStereoChemDoneProp)
{
    // SKETCH-1701

    TestSketcherWidget sk;
    sk.addFromString("C[C@H](N)O", Format::SMILES);
    auto molblock = sk.getString(Format::MDL_MOLV3000);
    BOOST_TEST(molblock.find(">  <_StereochemDone>") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(testPreserveStereoGroupIds)
{
    // SKETCH-2000

    TestSketcherWidget sk;
    sk.addFromString("C[C@H](O)Cl |o5:1|", Format::EXTENDED_SMILES);

    auto cxsmiles = sk.getString(Format::EXTENDED_SMILES);
    BOOST_CHECK(cxsmiles.find(" |o5:") != std::string::npos);

    auto molblock = sk.getString(Format::MDL_MOLV3000);
    BOOST_CHECK(molblock.find("M  V30 MDLV30/STEREL5 ATOMS=(1 ") !=
                std::string::npos);
}

BOOST_DATA_TEST_CASE(TestReadingMaeInputsWithInvalidChiralityLabels,
                     boost::unit_test::data::make(std::vector<std::string>{
                         "12_R_1_3_4_5",    // missing chiral atom
                         "12_ANR_1_3_4_15", // missing substituent atom
                         "2_S_1_3_4",       // incomplete substituent list
                         "2_ANS_1_3_4_5_2", // self bond and too many atoms
                     }),
                     invalid_chirality_label)
{
    constexpr const char* maeblock_template = R"({{
 s_m_m2io_version
 :::
 2.0.0
}}

f_m_ct {{
  i_m_ct_stereo_status
  s_m_title
  s_st_Chirality_1
  :::
  1
  ""
  {}
  m_atom[5] {{
    # First column is Index #
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord
    i_m_atomic_number
    i_m_formal_charge
    :::
    1 -1.000000 -0.060606 0.000000 6 0
    2 -1.000000 -1.560606 0.000000 6 0
    3 -0.250000 -2.859644 0.000000 17 0
    4 -2.500000 -1.560606 0.000000 9 0
    5 0.299038 -0.810606 0.000000 35 0
    :::
  }}
  m_bond[4] {{
    # First column is Index #
    i_m_from
    i_m_order
    i_m_to
    :::
    1 1 1 2
    2 2 1 3
    3 2 1 4
    4 2 1 5
    :::
  }}
}}
)";
    auto maeblock = fmt::format(maeblock_template, invalid_chirality_label);

    TestSketcherWidget sk;
    // we should be able to import the text without any issues
    sk.addFromString(maeblock, Format::MAESTRO);
    BOOST_TEST(sk.getString(Format::SMILES) == "CC(F)(Cl)Br");
}

BOOST_AUTO_TEST_CASE(test_wasm_API)
{
    TestSketcherWidget sk;

    auto assert_roundtrip = [&sk](auto rdkit_read, const auto& filename) {
        auto text = read_testfile(filename);
        // Body of sketcherImportText
        sk.clear();
        sk.addFromString(text, Format::AUTO_DETECT);
        // Body of sketcherExportMolBlock
        auto molblock = sk.getString(Format::MDL_MOLV3000);
        // Confirm RDKit can deserialize into a valid object
        auto rdkit_ptr = rdkit_read(molblock);
        BOOST_TEST(rdkit_ptr != nullptr);
    };

    auto read_mol = [](const auto& molblock) {
        return std::shared_ptr<RDKit::RWMol>(RDKit::MolBlockToMol(molblock));
    };

    // roundtrip a structure through the sketcher
    assert_roundtrip(read_mol, "structure_1.sdf");

    auto read_rxn = [](const auto& molblock) {
        return std::shared_ptr<RDKit::ChemicalReaction>(
            RDKit::RxnBlockToChemicalReaction(molblock));
    };

    // roundtrip a rection through the sketcher
    assert_roundtrip(read_rxn, "SKETCH_1003.rxn");
    assert_roundtrip(read_rxn, "rxn_no_r_group.rxn");
    assert_roundtrip(read_rxn, "rxn_r_group.rxn");
    assert_roundtrip(read_rxn, "rxn_aryl.rxn");
    assert_roundtrip(read_rxn, "reaction_smarts.rsmi");
}
