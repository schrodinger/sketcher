#define BOOST_TEST_MODULE Test_Sketcher

#include <memory>
#include <tuple>
#include <vector>

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <QGraphicsItem>
#include <QLineF>
#include <QPointF>
#include <QPolygonF>
#include <QtGlobal>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
// Boost doesn't know how to print QPoints or QStrings
BOOST_TEST_DONT_PRINT_LOG_VALUE(QPointF);
BOOST_TEST_DONT_PRINT_LOG_VALUE(QString);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

class TestBondItem : public BondItem
{
  public:
    using BondItem::calcArrowTip;
    using BondItem::calcDashedWedge;
    using BondItem::calcDashedWedgeParams;
    using BondItem::calcDoubleBondLines;
    using BondItem::calcSolidWedge;
    using BondItem::calcSquigglyWedge;
    using BondItem::calculateLinesToPaint;
    using BondItem::calcWedgeEnds;
    using BondItem::findBestRingForBond;
    using BondItem::findRingCenter;
    using BondItem::m_annotation_text;
    using BondItem::m_bond;

    TestBondItem(RDKit::Bond* bond, const AtomItem& start_atom_item,
                 const AtomItem& end_atom_item, Fonts& fonts,
                 const BondDisplaySettings& settings,
                 QGraphicsItem* parent = nullptr) :
        BondItem(bond, start_atom_item, end_atom_item, fonts, settings, parent)
    {
    }
};

std::pair<std::vector<std::shared_ptr<TestBondItem>>,
          std::shared_ptr<TestScene>>
createStructure(std::string input_text)
{
    auto test_scene = TestScene::getScene();
    import_mol_text(test_scene->m_mol_model, input_text);
    std::vector<AtomItem*> atom_items;
    std::vector<BondItem*> scene_bond_items;
    for (auto* item : test_scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_items.push_back(atom_item);
        }
        if (auto* bond_item = qgraphicsitem_cast<BondItem*>(item)) {
            scene_bond_items.push_back(bond_item);
        }
    }
    int idx = 0;
    std::vector<std::shared_ptr<TestBondItem>> bond_items;
    for (auto bond : test_scene->m_mol_model->getMol()->bonds()) {
        BOOST_TEST_REQUIRE(bond != nullptr);
        auto start_atom_idx = bond->getBeginAtomIdx();
        auto end_atom_idx = bond->getEndAtomIdx();
        auto bond_item = std::make_shared<TestBondItem>(
            bond, *atom_items[start_atom_idx], *atom_items[end_atom_idx],
            test_scene->m_fonts,
            *test_scene->m_sketcher_model->getBondDisplaySettingsPtr());
        bond_item->setPos(scene_bond_items.at(idx++)->pos());
        bond_items.push_back(bond_item);
    }
    return {bond_items, test_scene};
}

std::tuple<std::shared_ptr<TestBondItem>, std::shared_ptr<TestScene>, QLineF>
createBondItem()
{
    auto test_scene = TestScene::getScene();
    import_mol_text(test_scene->m_mol_model, "CC");
    std::vector<AtomItem*> atom_items;
    for (auto* item : test_scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_items.push_back(atom_item);
        }
    }
    BOOST_TEST_REQUIRE(atom_items.size() == 2);
    RDKit::Bond* bond = *(test_scene->m_mol_model->getMol()->bonds().begin());
    BOOST_TEST_REQUIRE(bond != nullptr);
    auto bond_item = std::make_shared<TestBondItem>(
        bond, *atom_items[0], *atom_items[1], test_scene->m_fonts,
        *test_scene->m_sketcher_model->getBondDisplaySettingsPtr());
    QLineF trimmed_line(10, 0, 10, 100);
    return std::make_tuple(bond_item, test_scene, trimmed_line);
}

BOOST_AUTO_TEST_CASE(test_calcArrowTip)
{
    auto [bond_item, test_scene, trimmed_line] = createBondItem();
    QPolygonF arrow_tip = bond_item->calcArrowTip(trimmed_line);
    BOOST_TEST(arrow_tip.size() == 4);
    BOOST_TEST(arrow_tip.isClosed());
}

BOOST_AUTO_TEST_CASE(test_chiral_bond_methods)
{
    auto [bond_item, test_scene, trimmed_line] = createBondItem();

    // test calcWedgeEnds
    {
        auto [wedge1, wedge2, wedge3] = bond_item->calcWedgeEnds(trimmed_line);
        BOOST_TEST(!wedge1.isNull());
        BOOST_TEST(!wedge2.isNull());
        BOOST_TEST(!wedge3.isNull());
        BOOST_TEST(wedge1 != wedge2);
        BOOST_TEST(wedge1 != wedge3);
        BOOST_TEST(wedge2 != wedge3);
    }

    // test calcSolidWedge
    {
        QPolygonF solid_wedge = bond_item->calcSolidWedge(trimmed_line);
        BOOST_TEST(solid_wedge.size() == 4);
        BOOST_TEST(solid_wedge.isClosed());
    }

    // test calcDashedWedge
    {
        auto dashes = bond_item->calcDashedWedge(trimmed_line);
        // dash.size() here is intentionally one greater than dash.size() from
        // calcSquigglyWedge.  See the comment in the calcDashedWedge
        // implementation for more details.
        BOOST_TEST(dashes.size() == 21);
        qreal length_of_previous_dash = -1;
        qreal y_of_previous_dash = -1;
        for (auto dash : dashes) {
            BOOST_TEST(dash.length() > length_of_previous_dash);
            BOOST_CHECK_CLOSE(dash.y1(), dash.y2(), 0.001);
            BOOST_TEST(dash.y1() > y_of_previous_dash);
            length_of_previous_dash = dash.length();
            y_of_previous_dash = dash.y1();
        }
    }

    // test calcSquigglyWedge
    {
        auto dashes = bond_item->calcSquigglyWedge(trimmed_line);
        BOOST_TEST(dashes.size() == 20);
        qreal length_of_previous_dash = -1;
        qreal y1_of_previous_dash = -1;
        qreal y2_of_previous_dash = -1;
        for (auto dash : dashes) {
            BOOST_TEST(dash.length() > length_of_previous_dash);
            BOOST_TEST(dash.y1() != dash.y2());
            bool one_side_further = dash.y1() > y1_of_previous_dash ||
                                    dash.y2() > y2_of_previous_dash;
            BOOST_TEST(one_side_further);
            length_of_previous_dash = dash.length();
            y1_of_previous_dash = dash.y1();
            y2_of_previous_dash = dash.y2();
        }
    }
}

BOOST_AUTO_TEST_CASE(test_findRingCenter)
{
    auto [bond_items, test_scene] = createStructure("C1CCCCC1");
    auto first_bond = bond_items[0];
    const RDKit::ROMol& molecule = first_bond->m_bond->getOwningMol();
    RDKit::RingInfo* ring_info = molecule.getRingInfo();
    auto ring_center = first_bond->mapToScene(
        first_bond->findRingCenter(molecule, ring_info, 0));

    // test that each bond of the ring finds the same ring center
    for (auto bond : bond_items) {
        auto ring_center_relative_to_bond =
            bond->findRingCenter(molecule, ring_info, 0);
        auto absolute_center = bond->mapToScene(ring_center_relative_to_bond);
        const auto SMALL = 0.01;
        BOOST_TEST((ring_center - absolute_center).manhattanLength() < SMALL);
    }
}

BOOST_AUTO_TEST_CASE(test_findBestRingForBond)
{
    auto [bond_items, test_scene] =
        createStructure("O=C(O)C1=NN=C(C2C3=CC=C(Cl)C=C3CCC3=C2N=CC=C3)C=C1");
    const RDKit::ROMol& molecule = bond_items[0]->m_bond->getOwningMol();
    RDKit::RingInfo* ring_info = molecule.getRingInfo();
    auto bond_rings = ring_info->bondRings();
    for (auto bond_item : bond_items) {
        auto bond = bond_item->m_bond;
        if (bond->getBondType() != RDKit::Bond::BondType::DOUBLE &&
            bond->getBondType() != RDKit::Bond::BondType::AROMATIC) {
            continue;
        }
        auto ring_idx = bond_item->findBestRingForBond(molecule, ring_info);
        // ignore the double bond in the carboxyl group, which isn't in any
        // rings
        if (ring_idx > 0) {
            // make sure that the bonds are drawn inside the 6 membered rings,
            // not the 7 membered ring
            auto ring_size = bond_rings[ring_idx].size();
            BOOST_TEST(ring_size == 6);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_stereo_label)
{
    // from a given input text, create sketcher objects and check that the
    // bond types, directions, and labels are correctly set for each bond
    auto check_bond_items = [](auto input_text, auto expected) {
        auto [bond_items, scene] = createStructure(input_text);
        BOOST_REQUIRE(bond_items.size() == expected.size());
        for (size_t i = 0; i < bond_items.size(); ++i) {
            auto bond_item = bond_items[i];
            auto bond = bond_item->m_bond;
            auto [bond_type, bond_dir, bond_annotation] = expected[i];
            BOOST_TEST(get_bond_type_and_query_label(bond).first == bond_type);
            BOOST_TEST(bond->getBondDir() == bond_dir);
            BOOST_TEST(bond_item->m_annotation_text.toStdString() ==
                       bond_annotation);
        }
    };

    using RDKit::Bond;
    std::vector<std::tuple<Bond::BondType, Bond::BondDir, std::string>>
        expected;

    // test that the stereo labels are correctly set for each bond of two
    // structures (non stereo bonds should have an empty label)
    expected = {{Bond::BondType::SINGLE, Bond::BondDir::ENDUPRIGHT, ""},
                {Bond::BondType::DOUBLE, Bond::BondDir::NONE, "(E)"},
                {Bond::BondType::SINGLE, Bond::BondDir::ENDUPRIGHT, ""}};
    check_bond_items("C\\C=C\\C", expected);

    expected = {{Bond::BondType::SINGLE, Bond::BondDir::ENDUPRIGHT, ""},
                {Bond::BondType::DOUBLE, Bond::BondDir::NONE, "(Z)"},
                {Bond::BondType::SINGLE, Bond::BondDir::ENDDOWNRIGHT, ""}};
    check_bond_items("C\\C=C/C", expected);

    // test that stereo wedging is correctly set on each bond; for input SMILES
    // wedging is explicitly calculated by the sketcher code from the input
    // parities detected by RDKit on SMILES read
    expected = {{Bond::BondType::SINGLE, Bond::BondDir::NONE, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::NONE, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::NONE, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::BEGINDASH, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::NONE, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::BEGINWEDGE, ""},
                {Bond::BondType::SINGLE, Bond::BondDir::NONE, ""}};
    check_bond_items("CC(C)[C@@H](C)[C@H](C)N", expected);

    // test the above with MDL input where the bond directions are explicitly
    // set in the mol block (as opposed to calculated from parities)
    std::swap(expected[3], expected[5]); // the wedge/dash are swapped in MDL!
    std::string molblock{R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.301429 6.542858 0.000000 0
M  V30 2 C -1.304286 5.114286 0.000000 0
M  V30 3 C -2.542857 4.402572 0.000000 0
M  V30 4 C -0.068572 4.397429 0.000000 0
M  V30 5 C 1.170000 5.109429 0.000000 0
M  V30 6 C -0.071429 2.968857 0.000000 0
M  V30 7 C 1.164285 2.252286 0.000000 0
M  V30 8 N -1.310286 2.257143 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 4 5 CFG=1
M  V30 5 1 4 6
M  V30 6 1 6 7 CFG=3
M  V30 7 1 6 8
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 4 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
)MDL"};
    check_bond_items(molblock, expected);
}

/**
 * test that the dashed wedge is calculated for a bond with a zero length
 * without crashing and that the result is a single dash (SKETCH-2148)
 */
BOOST_AUTO_TEST_CASE(test_dashed_bond_zero_length)
{
    auto [bond_item, test_scene, trimmed_line] = createBondItem();
    auto [num_dashes, wedge_start, offset_towards_p2, offset_towards_p3,
          dashes] = bond_item->calcDashedWedgeParams(QLineF(0, 0, 0, 0));
    BOOST_TEST(num_dashes == 1);
}

BOOST_AUTO_TEST_CASE(test_bond_stereo_tooltips)
{
    // Test that E/Z double bonds have correct tooltips
    // The annotation text comes with parentheses like "(E)" or "(Z)"
    // and the tooltip should preserve them for consistency with atom stereo
    // labels
    auto [bond_items, scene] = createStructure("C/C(N)=C(/C)O");

    // Find the double bond with E stereochemistry
    bool found_stereo_bond = false;
    for (auto bond_item : bond_items) {
        if (bond_item->m_annotation_text == "(E)") {
            BOOST_TEST(bond_item->toolTip() == "Stereo: (E)");
            found_stereo_bond = true;
        }
    }
    BOOST_TEST(found_stereo_bond);

    // Test Z stereochemistry
    auto [bond_items_z, scene_z] = createStructure("C\\C=C/C");
    found_stereo_bond = false;
    for (auto bond_item : bond_items_z) {
        if (bond_item->m_annotation_text == "(Z)") {
            BOOST_TEST(bond_item->toolTip() == "Stereo: (Z)");
            found_stereo_bond = true;
        }
    }
    BOOST_TEST(found_stereo_bond);
}

} // namespace sketcher
} // namespace schrodinger
