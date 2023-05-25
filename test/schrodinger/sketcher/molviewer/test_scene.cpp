#define BOOST_TEST_MODULE Test_Sketcher

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
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
    TestScene test_scene;
    BOOST_TEST(test_scene.fontSize() == DEFAULT_FONT_SIZE);
    test_scene.setFontSize(8);
    BOOST_TEST(test_scene.fontSize() == 8);
    test_scene.setFontSize(45);
    BOOST_TEST(test_scene.fontSize() == 45);
}

void count_visible_atoms(const Scene& test_scene, unsigned& num_visible_atoms,
                         unsigned& num_hidden_atoms)
{
    num_visible_atoms = num_hidden_atoms = 0;
    for (auto item : test_scene.items()) {
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
    TestScene test_scene;
    unsigned num_visible_atoms = 0;
    unsigned num_hidden_atoms = 0;
    test_scene.m_mol_model->addMolFromText("CCCCC", Format::SMILES);

    // all carbons should be hidden
    test_scene.setCarbonsLabeled(CarbonLabels::NONE);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 0);
    BOOST_TEST(num_hidden_atoms == 5);

    // only terminal carbons should be visible
    test_scene.setCarbonsLabeled(CarbonLabels::TERMINAL);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 2);
    BOOST_TEST(num_hidden_atoms == 3);

    // all carbons should be visible
    test_scene.setCarbonsLabeled(CarbonLabels::ALL);
    count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms);
    BOOST_TEST(num_visible_atoms == 5);
    BOOST_TEST(num_hidden_atoms == 0);
}

BOOST_AUTO_TEST_CASE(test_getInteractiveItems)
{
    TestScene scene;
    auto items = scene.getInteractiveItems();
    BOOST_TEST(items.size() == 0);
    BOOST_TEST(scene.items().size() == 3); // selection items

    scene.m_mol_model->addMolFromText("CC", Format::SMILES);
    items = scene.getInteractiveItems();
    BOOST_TEST(items.size() == 3); // two atoms and a bond
    BOOST_TEST(scene.items().size() == 6);

    scene.m_mol_model->clear();
    items = scene.getInteractiveItems();
    BOOST_TEST(items.size() == 0);
    BOOST_TEST(scene.items().size() == 3);
}

BOOST_AUTO_TEST_CASE(test_item_selection)
{
    TestScene scene;
    scene.m_mol_model->addMolFromText("CC", Format::SMILES);

    auto test_selected_items = [](const auto& scene, bool expected) {
        BOOST_TEST(scene.m_selection_highlighting_item->isVisible() ==
                   expected);
        for (auto item : scene.getInteractiveItems()) {
            BOOST_TEST(item->flags() & QGraphicsItem::ItemIsSelectable);
            BOOST_TEST(item->isSelected() == expected);
        }
    };

    // no selection
    test_selected_items(scene, false);

    // select everything
    scene.m_mol_model->selectAll();
    test_selected_items(scene, true);

    // clear selection
    scene.m_mol_model->clearSelection();
    test_selected_items(scene, false);

    // invert selection
    scene.m_mol_model->invertSelection();
    test_selected_items(scene, true);
    scene.m_mol_model->invertSelection();
    test_selected_items(scene, false);
}

} // namespace sketcher
} // namespace schrodinger
