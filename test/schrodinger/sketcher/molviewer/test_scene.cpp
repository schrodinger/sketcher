#define BOOST_TEST_MODULE Test_Sketcher

#include "schrodinger/sketcher/molviewer/scene.h"

#include <QRectF>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "schrodinger/rdkit_extensions/convert.h"

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_importText)
{
    Scene test_scene;
    test_scene.importText("c1nccc2n1ccc2", Format::SMILES);
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

    test_scene.clear();
    mol = test_scene.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 0);

    // import failed, exception caught, still an empty scene
    test_scene.importText("nonsense", Format::AUTO_DETECT);
    mol = test_scene.getRDKitMolecule();
    BOOST_TEST(mol->getNumAtoms() == 0);
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
    Scene test_scene;
    unsigned num_visible_atoms = 0;
    unsigned num_hidden_atoms = 0;
    test_scene.importText("CCCCC", Format::SMILES);

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

} // namespace sketcher
} // namespace schrodinger
