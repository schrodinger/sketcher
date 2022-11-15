#define BOOST_TEST_MODULE Test_Sketcher
#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <QRectF>

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_load_smiles)
{
    Scene test_scene;
    test_scene.loadSmiles("c1nccc2n1ccc2");
    auto mol = test_scene.getRDKitMolecule();
    BOOST_TEST_REQUIRE(mol != nullptr);
    BOOST_TEST(mol->getNumAtoms() == 9);
    unsigned num_atoms = 0;
    unsigned num_bonds = 0;
    for (auto item : test_scene.items()) {
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
}

BOOST_AUTO_TEST_CASE(test_load_mol)
{
    Scene test_scene;
    std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCCC"));
    BOOST_TEST_REQUIRE(mol != nullptr);
    RDDepict::compute2DCoords(*mol, nullptr, true);
    // pass mol to loadMol
    test_scene.loadMol(*mol);
    BOOST_TEST(test_scene.getRDKitMolecule()->getNumAtoms() == 4);
    // pass shared_ptr to loadMol
    test_scene.loadMol(mol);
    BOOST_TEST(test_scene.getRDKitMolecule()->getNumAtoms() == 4);
}

BOOST_AUTO_TEST_CASE(test_font_size)
{
    Scene test_scene;
    BOOST_TEST(test_scene.fontSize() == DEFAULT_FONT_SIZE);
    test_scene.setFontSize(8);
    BOOST_TEST(test_scene.fontSize() == 8);
    test_scene.setFontSize(45);
    BOOST_TEST(test_scene.fontSize() == 45);
}

bool count_visible_atoms(const Scene& test_scene, unsigned& num_visible_atoms,
                         unsigned& num_hidden_atoms)
{
    num_visible_atoms = num_hidden_atoms = 0;
    for (auto item : test_scene.items()) {
        if (item->type() == AtomItem::Type) {
            QRectF brect = item->boundingRect();
            if (brect.width() == 0 && brect.height() == 0) {
                ++num_hidden_atoms;
            } else if (brect.width() > 0 && brect.height() > 0) {
                ++num_visible_atoms;
            } else {
                // The bounding rect should never be a line
                return false;
            }
        }
    }
    return true;
}

BOOST_AUTO_TEST_CASE(test_all_atoms_shown)
{
    Scene test_scene;
    unsigned num_visible_atoms = 0;
    unsigned num_hidden_atoms = 0;
    test_scene.loadSmiles("CCCCC");

    // all carbons should be hidden
    test_scene.setCarbonsLabeled(CarbonLabels::NONE);
    BOOST_TEST(
        count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms));
    BOOST_TEST(num_visible_atoms == 0);
    BOOST_TEST(num_hidden_atoms == 5);

    // only terminal carbons should be visible
    test_scene.setCarbonsLabeled(CarbonLabels::TERMINAL);
    BOOST_TEST(
        count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms));
    BOOST_TEST(num_visible_atoms == 2);
    BOOST_TEST(num_hidden_atoms == 3);

    // all carbons should be visible
    test_scene.setCarbonsLabeled(CarbonLabels::ALL);
    BOOST_TEST(
        count_visible_atoms(test_scene, num_visible_atoms, num_hidden_atoms));
    BOOST_TEST(num_visible_atoms == 5);
    BOOST_TEST(num_hidden_atoms == 0);
}

} // namespace sketcher
} // namespace schrodinger
