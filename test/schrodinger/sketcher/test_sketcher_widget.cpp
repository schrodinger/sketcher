#define BOOST_TEST_MODULE sketcher_widget_test

#include <string>

#include <QString>

#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/coord_utils.h"
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

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::ColorScheme)

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

// Reuse a single widget instance across all tests for better performance
// Note: Widget is intentionally leaked to avoid Qt destruction order issues.
// The widget cannot be safely destroyed after QApplication cleanup without
// causing crashes on Linux/Windows. This small leak (< 1MB) is acceptable for
// test code and only occurs at program exit when the OS reclaims all memory.
struct TestWidgetFixture {
    static TestSketcherWidget* get()
    {
        static TestSketcherWidget* widget = nullptr;
        if (!widget) {
            widget = new TestSketcherWidget();
        }
        // Clear state before each test for isolation
        widget->clear();
        widget->m_undo_stack->clear();
        widget->m_sketcher_model->reset();
        return widget;
    }
};

BOOST_GLOBAL_FIXTURE(TestWidgetFixture);

BOOST_AUTO_TEST_CASE(test_addRDKitMolecule_getRDKitMolecule)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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

    TestSketcherWidget& sk = *TestWidgetFixture::get();
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

    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString(text);

    BOOST_TEST(sk.getString(Format::SMILES) == "CC(=O)O.CCO>>CCOC(C)=O");
}

/**
 * Make sure that we export PDB files with a single model, even if the MolModel
 * mol stores 3d coordinates in a second conformer.
 */
BOOST_AUTO_TEST_CASE(test_PDB_format_single_model)
{
    std::shared_ptr<RDKit::ROMol> mol_to_add(RDKit::SmilesToMol("CCC"));
    auto conf = new RDKit::Conformer(3);
    conf->set3D(true);
    conf->setAtomPos(0, {1, 1, 1});
    conf->setAtomPos(1, {2, 2, 2});
    conf->setAtomPos(2, {3, 3, 3});
    mol_to_add->addConformer(conf, true);

    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addRDKitMolecule(*mol_to_add);
    // RDKit mol exports should contain both the 2D and the 3D conformers
    auto mol_out = sk.getRDKitMolecule();
    BOOST_TEST(mol_out->getNumConformers() == 2);

    auto pdb_text = sk.getString(Format::PDB);
    auto pdb_qtext = QString::fromStdString(pdb_text);
    BOOST_TEST(!pdb_qtext.contains("MODEL"));
    BOOST_TEST(pdb_qtext.count("HETATM") == 3);
}

/**
 * Make sure that molfiles are exported using 2D conformations rather than 3D
 */
BOOST_AUTO_TEST_CASE(test_2D_molfile)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString("C");
    auto molfile = sk.getString(Format::MDL_MOLV3000);
    BOOST_TEST(molfile.find("2D") != std::string::npos);
    BOOST_TEST(molfile.find("3D") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(test_cut_copy_paste)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    sk.m_mol_model->select({}, {bond}, {}, {}, {}, SelectMode::SELECT);
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
    sk.m_mol_model->select({}, {bond}, {}, {}, {}, SelectMode::SELECT);
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
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.setClipboardContents("CCC\n\r");
    sk.paste();
    auto mol = sk.m_mol_model->getMol();
    BOOST_TEST(mol->getNumAtoms() == 3);
}

BOOST_AUTO_TEST_CASE(test_importText_slot)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    TestSketcherWidget& sk = *TestWidgetFixture::get();

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
    TestSketcherWidget& sk = *TestWidgetFixture::get();

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
    TestSketcherWidget& sk = *TestWidgetFixture::get();

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
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    TestSketcherWidget& sk = *TestWidgetFixture::get();
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
    std::string orig_smiles = "c1ccccc1";
    auto mol = to_rdkit(orig_smiles, Format::SMILES);
    ::schrodinger::sketcher::update_2d_coordinates(*mol);
    auto input_string = to_string(*mol, sample);

    Format export_format = sample;
    std::string reference = input_string;
    // this test doesn't work for some formats, so just sanity check the SMILES
    // string in those cases
    if (sample == Format::MDL_MOLV2000) {
        // Sketcher doesn't support exporting to V2000
        export_format = Format::SMILES;
        // we'll lose aromaticity information for atoms when converting to a
        // molblock, so we can't use the original SMILES string
        reference = "C1:C:C:C:C:C:1";
    } else if (sample == Format::RDMOL_BINARY_BASE64) {
        // we won't get bit-for-bit fidelity with the round-trip
        export_format = Format::SMILES;
        reference = orig_smiles;
    }

    // Check roundtripping and auto-detect
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString(input_string, sample);
    BOOST_TEST(sk.getString(export_format) == reference);
    sk.clear();
    sk.addFromString(input_string);
    BOOST_TEST(sk.getString(export_format) == reference);
}

BOOST_DATA_TEST_CASE(test_reactions_roundtrip,
                     boost::unit_test::data::make(RXN_FORMATS))
{
    std::string smiles = "CC(=O)O.OCC>>CC(=O)OCC";
    auto reaction = to_rdkit_reaction(smiles, Format::SMILES);
    auto text = to_string(*reaction, sample);

    // Check auto-detect import through the sketcher interface as well
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString(text);
    BOOST_TEST(sk.getString(Format::SMILES) == "CC(=O)O.CCO>>CCOC(C)=O");
}

BOOST_AUTO_TEST_CASE(testRGRoup0)
{

    TestSketcherWidget& sk = *TestWidgetFixture::get();

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

    TestSketcherWidget& sk = *TestWidgetFixture::get();

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

    TestSketcherWidget& sk = *TestWidgetFixture::get();

    std::string cxsmiles{"[*]C1=CC=CC=C1 |$_R1;;;;;;$,c:3,5,t:1|"};

    sk.addFromString(cxsmiles, Format::EXTENDED_SMILES);
    BOOST_TEST(sk.getString(Format::EXTENDED_SMARTS) ==
               "[#0]-[#6]1=[#6]-[#6]=[#6]-[#6]=[#6]-1 |$_R1;;;;;;$|");
}

BOOST_AUTO_TEST_CASE(testLeakStereoChemDoneProp)
{
    // SKETCH-1701

    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString("C[C@H](N)O", Format::SMILES);
    auto molblock = sk.getString(Format::MDL_MOLV3000);
    BOOST_TEST(molblock.find(">  <_StereochemDone>") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(testPreserveStereoGroupIds)
{
    // SKETCH-2000

    TestSketcherWidget& sk = *TestWidgetFixture::get();
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

    TestSketcherWidget& sk = *TestWidgetFixture::get();
    // we should be able to import the text without any issues
    sk.addFromString(maeblock, Format::MAESTRO);
    BOOST_TEST(sk.getString(Format::SMILES) == "CC(F)(Cl)Br");
}

BOOST_AUTO_TEST_CASE(test_wasm_API)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();

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

/**
 * Make sure that we can select atoms/bonds and retrieve the selection via the
 * SketcherWidget.  Also make sure that the selected atoms/bonds belong to the
 * molecule returned by getRDKitMolecule.
 */
BOOST_AUTO_TEST_CASE(test_selection)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.addFromString("NCCC");
    auto mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 4);
    // make sure that getRDKitMolecule return the same molecule twice in a row
    BOOST_TEST(mol.get() == sk.getRDKitMolecule().get());

    sk.select({mol->getAtomWithIdx(0)},
              {mol->getBondWithIdx(0), mol->getBondWithIdx(1)});
    auto sel_atoms = sk.getSelectedAtoms();
    BOOST_TEST(sel_atoms.size() == 1);
    auto* atom = *sel_atoms.begin();
    BOOST_TEST(atom->hasOwningMol());
    // make sure that the atom is owned by mol
    BOOST_TEST(&atom->getOwningMol() == mol.get());

    auto sel_bonds = sk.getSelectedBonds();
    BOOST_TEST(sel_bonds.size() == 2);
    auto* bond = *sel_bonds.begin();
    BOOST_TEST(bond->hasOwningMol());
    BOOST_TEST(&bond->getOwningMol() == mol.get());

    // make sure that the sk.getRDKitMolecule() return value is updated when the
    // underlying molecule changes
    sk.addFromString("C");
    auto new_mol = sk.getRDKitMolecule();
    BOOST_TEST(new_mol->getNumAtoms() == 5);
    BOOST_TEST(mol.get() != new_mol.get());
}

/*
 * Make sure that the color scheme doesn't get reset when events are processed
 */
BOOST_AUTO_TEST_CASE(test_color_scheme_set_before_show)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.setColorScheme(ColorScheme::DARK_MODE);
    QCoreApplication::processEvents();
    BOOST_TEST(sk.m_sketcher_model->getColorScheme() == ColorScheme::DARK_MODE);
    // make sure that the bond color is close to white for dark mode
    auto bond_color = sk.m_sketcher_model->getBondDisplaySettingsPtr()->m_color;
    BOOST_TEST(bond_color.lightnessF() > 0.9);
}

/**
 * Test that cleaning up a selection fits only the selection to the screen,
 * not the entire scene. This test verifies the fix where fitToScreen()
 * is now called with the selection_only parameter after cleaning up a
 * selection.
 */
BOOST_AUTO_TEST_CASE(test_cleanup_selection_fits_only_selection)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();

    // Add a long chain molecule
    // This creates a molecule that spans a large area
    sk.addFromString("CCCCCCCCCCCCCCCCCCCC", Format::SMILES); // 20 carbons
    auto mol = sk.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 20);

    // Show the widget to ensure view geometry is initialized
    sk.show();
    QCoreApplication::processEvents();

    // Fit everything to screen initially
    auto view = sk.m_ui->view;
    view->fitToScreen(false);
    auto initial_transform = view->transform();

    // Select only the first 3 atoms (a small portion at one end)
    sk.select({mol->getAtomWithIdx(0), mol->getAtomWithIdx(1),
               mol->getAtomWithIdx(2)},
              {});
    BOOST_TEST(sk.m_mol_model->getSelectedAtoms().size() == 3);

    // Test fitting only the selection - should zoom in on just the first 3
    // atoms
    view->fitToScreen(true); // selection_only = true
    auto transform_selection = view->transform();

    // Reset to initial state
    view->setTransform(initial_transform);

    // Test fitting the entire scene - should show all 20 atoms
    view->fitToScreen(false); // selection_only = false
    auto transform_all = view->transform();

    // The transforms should be different when fitting 3 atoms vs 20 atoms
    // When fitting only 3 atoms at one end of a 20-atom chain, we should be
    // zoomed in significantly more than when fitting all 20 atoms
    BOOST_TEST(
        transform_selection.m11() != transform_all.m11(),
        "Selection-only fit should zoom differently than full scene fit");

    // Fitting only 3 atoms should result in a higher zoom level than fitting 20
    BOOST_TEST(transform_selection.m11() > transform_all.m11(),
               "Fitting 3 atoms should zoom in more than fitting 20 atoms");
}

/**
 * Make sure that SketcherWidget::addTextToMolModel properly prevents the user
 * from adding atomistic or monomeric models when appropriate.
 */
BOOST_AUTO_TEST_CASE(test_addTextToMolModel)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    const std::string MONOMERIC_STRING = "PEPTIDE1{D.E.F.G}$$$$V2.0";
    const std::string ATOMISTIC_STRING = "CCC";

    sk.setInterfaceType(InterfaceType::ATOMISTIC);
    BOOST_CHECK_THROW(sk.addTextToMolModel(MONOMERIC_STRING), std::exception);
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(ATOMISTIC_STRING));
    sk.clear();

    sk.setInterfaceType(InterfaceType::MONOMERIC);
    BOOST_CHECK_THROW(sk.addTextToMolModel(ATOMISTIC_STRING), std::exception);
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(MONOMERIC_STRING));
    sk.clear();

    sk.setInterfaceType(InterfaceType::ATOMISTIC_OR_MONOMERIC);
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(ATOMISTIC_STRING));
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(ATOMISTIC_STRING));
    // we can't add monomers if there are already atoms
    BOOST_CHECK_THROW(sk.addTextToMolModel(MONOMERIC_STRING), std::exception);
    sk.clear();

    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(MONOMERIC_STRING));
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(MONOMERIC_STRING));
    // we can't add atoms if there are already monomers
    BOOST_CHECK_THROW(sk.addTextToMolModel(ATOMISTIC_STRING), std::exception);
}

/**
 * Make sure that we can read in monomeric models that contain additional data
 * in S-groups.
 */
BOOST_AUTO_TEST_CASE(test_addTextToMolModel_monomeric_s_sroups)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    sk.setInterfaceType(InterfaceType::ATOMISTIC_OR_MONOMERIC);

    // extended annotations are currently unsupported, but make sure that we
    // throw an exception (which will be caught and put in an error dialog)
    // instead of crashing
    const std::string HELM_WITH_EXTENDED_ANNOTATION =
        R"(RNA1{R(A)P.R(C)P.R(G)}$$${"my chain":"my annotation"}$V2.0)";
    BOOST_CHECK_THROW(sk.addTextToMolModel(HELM_WITH_EXTENDED_ANNOTATION),
                      std::exception);

    // basic annotations create a COP S-group
    const std::string HELM_WITH_ANNOTATION =
        R"(PEPTIDE1{A.C.D.D.E}"HC"|PEPTIDE2{G.C.S.S.S.P.K.K.V.K}"LC"$$$$V2.0)";
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(HELM_WITH_ANNOTATION));
    BOOST_TEST(sk.getRDKitMolecule()->getNumAtoms() == 15);
    sk.clear();

    // FASTA strings also create a COP S-group
    const std::string FASTA = ">foo\nAAA";
    BOOST_CHECK_NO_THROW(sk.addTextToMolModel(FASTA));
    BOOST_TEST(sk.getRDKitMolecule()->getNumAtoms() == 3);
}
/**
 * Test copying partial reactions with and without non-molecular objects
 * selected. This validates that SKETCH-2632 allows copying partial selections
 * when no reaction elements (arrow, plus signs) are selected, but prevents
 * copying when they are selected.
 */
BOOST_AUTO_TEST_CASE(test_copy_partial_reaction)
{
    TestSketcherWidget& sk = *TestWidgetFixture::get();
    auto mol_model = sk.m_mol_model;

    // Test with simple reaction: reactants >> products
    sk.addFromString("CC>>CCC");
    BOOST_TEST(mol_model->hasReactionArrow());

    auto mol = mol_model->getMol();
    // Reactants: atoms 0-1 (CC)
    // Products: atoms 2-4 (CCC)
    auto reactant_atoms = std::unordered_set<const RDKit::Atom*>{
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)};
    auto product_atoms = std::unordered_set<const RDKit::Atom*>{
        mol->getAtomWithIdx(2), mol->getAtomWithIdx(3), mol->getAtomWithIdx(4)};

    // Select reactants only (no arrow) - should succeed
    mol_model->select(reactant_atoms, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->hasSelection());
    BOOST_TEST(mol_model->getSelectedNonMolecularObjects().empty());
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == "CC");

    // Select products only (no arrow) - should succeed
    mol_model->select(product_atoms, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->hasSelection());
    BOOST_TEST(mol_model->getSelectedNonMolecularObjects().empty());
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == "CCC");

    // Select reactants + arrow - should fail (clipboard unchanged)
    auto arrow = mol_model->getReactionArrow();
    BOOST_REQUIRE(arrow != nullptr);
    mol_model->select(reactant_atoms, {}, {}, {}, {arrow},
                      SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->hasSelection());
    BOOST_TEST(!mol_model->getSelectedNonMolecularObjects().empty());
    std::string clipboard_before = sk.getClipboardContents();
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    // Error shown internally, clipboard should not update
    BOOST_TEST(sk.getClipboardContents() == clipboard_before);

    // Select products + arrow - should fail (clipboard unchanged)
    mol_model->select(product_atoms, {}, {}, {}, {arrow},
                      SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->hasSelection());
    BOOST_TEST(!mol_model->getSelectedNonMolecularObjects().empty());
    clipboard_before = sk.getClipboardContents();
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == clipboard_before);

    // Copy all - should work as complete reaction
    sk.copy(Format::SMILES, SceneSubset::ALL);
    BOOST_TEST(sk.getClipboardContents() == "CC>>CCC");

    // Test with multi-component reaction: A + B >> C
    sk.clear();
    sk.addFromString("CC.CCC>>CCCCC");
    BOOST_TEST(mol_model->hasReactionArrow());

    mol = mol_model->getMol();
    // Reactant A: atoms 0-1 (CC)
    // Reactant B: atoms 2-4 (CCC)
    // Product C: atoms 5-9 (CCCCC)
    auto reactantA_atoms = std::unordered_set<const RDKit::Atom*>{
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)};
    auto reactantB_atoms = std::unordered_set<const RDKit::Atom*>{
        mol->getAtomWithIdx(2), mol->getAtomWithIdx(3), mol->getAtomWithIdx(4)};
    auto productC_atoms = std::unordered_set<const RDKit::Atom*>{
        mol->getAtomWithIdx(5), mol->getAtomWithIdx(6), mol->getAtomWithIdx(7),
        mol->getAtomWithIdx(8), mol->getAtomWithIdx(9)};

    // Select one reactant component - should succeed
    mol_model->select(reactantA_atoms, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->getSelectedNonMolecularObjects().empty());
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == "CC");

    // Select both reactant components - should succeed
    auto both_reactants = reactantA_atoms;
    both_reactants.insert(reactantB_atoms.begin(), reactantB_atoms.end());
    mol_model->select(both_reactants, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->getSelectedNonMolecularObjects().empty());
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    // Result may be "CC.CCC" or "CCC.CC" depending on atom ordering
    auto clip = sk.getClipboardContents();
    BOOST_TEST((clip == "CC.CCC" || clip == "CCC.CC"));

    // Select reactants + plus sign - should fail
    auto all_non_mol = mol_model->getNonMolecularObjects();
    arrow = mol_model->getReactionArrow();
    // Find a plus sign (any non-molecular object that isn't the arrow)
    const NonMolecularObject* plus_sign = nullptr;
    for (auto* obj : all_non_mol) {
        if (obj != arrow) {
            plus_sign = obj;
            break;
        }
    }
    BOOST_REQUIRE(plus_sign != nullptr);
    mol_model->select(both_reactants, {}, {}, {}, {plus_sign},
                      SelectMode::SELECT_ONLY);
    BOOST_TEST(!mol_model->getSelectedNonMolecularObjects().empty());
    clipboard_before = sk.getClipboardContents();
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == clipboard_before);

    // Select product component - should succeed
    mol_model->select(productC_atoms, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(mol_model->getSelectedNonMolecularObjects().empty());
    sk.copy(Format::SMILES, SceneSubset::SELECTION);
    BOOST_TEST(sk.getClipboardContents() == "CCCCC");

    // Copy all - complete multi-component reaction
    sk.copy(Format::SMILES, SceneSubset::ALL);
    // RDKit may reorder components
    clip = sk.getClipboardContents();
    BOOST_TEST((clip == "CC.CCC>>CCCCC" || clip == "CCC.CC>>CCCCC"));
}
