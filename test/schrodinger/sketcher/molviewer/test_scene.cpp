#define BOOST_TEST_MODULE Test_Sketcher

#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <QRectF>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_font_size)
{
    auto test_scene = TestScene::getScene();
    BOOST_TEST(test_scene->m_fonts.size() == DEFAULT_FONT_SIZE);
    test_scene->setFontSize(8);
    BOOST_TEST(test_scene->m_fonts.size() == 8);
    test_scene->setFontSize(45);
    BOOST_TEST(test_scene->m_fonts.size() == 45);
}

void count_visible_atoms(const std::shared_ptr<Scene> test_scene,
                         unsigned& num_visible_atoms,
                         unsigned& num_hidden_atoms)
{
    num_visible_atoms = num_hidden_atoms = 0;
    for (auto item : test_scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            if (atom_item->labelIsVisible()) {
                ++num_visible_atoms;
            } else {
                ++num_hidden_atoms;
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_all_atoms_shown)
{
    auto test_scene = TestScene::getScene();
    unsigned num_visible_atoms = 0;
    unsigned num_hidden_atoms = 0;
    import_mol_text(test_scene->m_mol_model, "CCCCC", Format::SMILES);

    // all carbons should be hidden
    test_scene->setCarbonsLabeled(CarbonLabels::NONE);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 0);
    BOOST_TEST(num_hidden_atoms == 5);

    // only terminal carbons should be visible
    test_scene->setCarbonsLabeled(CarbonLabels::TERMINAL);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 2);
    BOOST_TEST(num_hidden_atoms == 3);

    // all carbons should be visible
    test_scene->setCarbonsLabeled(CarbonLabels::ALL);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 5);
    BOOST_TEST(num_hidden_atoms == 0);
}

BOOST_AUTO_TEST_CASE(test_getInteractiveItems)
{
    auto scene = TestScene::getScene();
    auto items = scene->getInteractiveItems();
    BOOST_TEST(items.size() == 0);
    BOOST_TEST(scene->items().size() == 4); // selection and highlighting items

    import_mol_text(scene->m_mol_model, "CC", Format::SMILES);
    items = scene->getInteractiveItems();
    BOOST_TEST(items.size() == 3); // two atoms and a bond
    BOOST_TEST(scene->items().size() == 7);

    scene->m_mol_model->clear();
    items = scene->getInteractiveItems();
    BOOST_TEST(items.size() == 0);
    BOOST_TEST(scene->items().size() == 4);
}

BOOST_AUTO_TEST_CASE(test_item_selection)
{
    auto scene = TestScene::getScene();
    import_mol_text(scene->m_mol_model, "CC", Format::SMILES);

    auto test_selected_items = [](const auto& scene, bool expected) {
        BOOST_TEST(scene->m_selection_highlighting_item->isVisible() ==
                   expected);
        for (auto item : scene->getInteractiveItems()) {
            BOOST_TEST(item->flags() & QGraphicsItem::ItemIsSelectable);
            BOOST_TEST(item->isSelected() == expected);
        }
    };

    // no selection
    test_selected_items(scene, false);

    // select everything
    scene->m_mol_model->selectAll();
    test_selected_items(scene, true);

    // clear selection
    scene->m_mol_model->clearSelection();
    test_selected_items(scene, false);

    // invert selection
    scene->m_mol_model->invertSelection();
    test_selected_items(scene, true);
    scene->m_mol_model->invertSelection();
    test_selected_items(scene, false);
}

BOOST_AUTO_TEST_CASE(test_ensureCompleteAttachmentPoints)
{
    auto scene = TestScene::getScene();
    import_mol_text(scene->m_mol_model, "CC* |$;;_AP1$|", Format::SMILES);
    auto mol = scene->m_mol_model->getMol();
    auto* ap_atom = mol->getAtomWithIdx(2);
    auto* ap_atom_item = scene->m_atom_to_atom_item.at(ap_atom);
    auto* ap_bond = mol->getBondWithIdx(1);
    auto* ap_bond_item = scene->m_bond_to_bond_item.at(ap_bond);

    // test with attachment point atom
    auto complete_items = scene->ensureCompleteAttachmentPoints({ap_atom_item});
    BOOST_TEST(complete_items.size() == 2);
    BOOST_TEST(complete_items.contains(ap_atom_item));
    BOOST_TEST(complete_items.contains(ap_bond_item));

    // test with attachment point bond
    complete_items = scene->ensureCompleteAttachmentPoints({ap_bond_item});
    BOOST_TEST(complete_items.size() == 2);
    BOOST_TEST(complete_items.contains(ap_atom_item));
    BOOST_TEST(complete_items.contains(ap_bond_item));

    // test with non-attachment point atom
    auto* c_atom = mol->getAtomWithIdx(1);
    auto* c_atom_item = scene->m_atom_to_atom_item.at(c_atom);
    complete_items = scene->ensureCompleteAttachmentPoints({c_atom_item});
    BOOST_TEST(complete_items.size() == 1);
    BOOST_TEST(complete_items.contains(c_atom_item));

    // test with non-attachment point bond
    auto* cc_bond = mol->getBondWithIdx(0);
    auto* cc_bond_item = scene->m_bond_to_bond_item.at(cc_bond);
    complete_items = scene->ensureCompleteAttachmentPoints({cc_bond_item});
    BOOST_TEST(complete_items.size() == 1);
    BOOST_TEST(complete_items.contains(cc_bond_item));
}

} // namespace sketcher
} // namespace schrodinger
