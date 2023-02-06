#define BOOST_TEST_MODULE Test_Sketcher

#include "schrodinger/sketcher/molviewer/bond_item.h"

#include <memory>
#include <tuple>
#include <vector>

#include <QtGlobal>
#include <QLineF>
#include <QGraphicsItem>
#include <QPolygonF>
#include <QPointF>

#include <GraphMol/CoordGen.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/scene.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);
// Boost doesn't know how to print QPoints
BOOST_TEST_DONT_PRINT_LOG_VALUE(QPointF);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

class TestScene : public Scene
{
  public:
    using Scene::m_bond_item_settings;
    using Scene::m_mol;
};

class TestBondItem : public BondItem
{
  public:
    using BondItem::arePointsOnSameSideOfLine;
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
    using BondItem::m_bond;

    TestBondItem(RDKit::Bond* bond, const AtomItem& start_atom_item,
                 const AtomItem& end_atom_item, BondItemSettings& settings,
                 QGraphicsItem* parent = nullptr) :
        BondItem(bond, start_atom_item, end_atom_item, settings, parent)
    {
    }
};

std::pair<std::vector<std::shared_ptr<TestBondItem>>,
          std::shared_ptr<TestScene>>
createStructure(std::string smiles)
{
    auto test_scene = std::make_shared<TestScene>();
    test_scene->importText(smiles, Format::SMILES);
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
    for (auto bond : test_scene->m_mol->bonds()) {
        BOOST_TEST_REQUIRE(bond != nullptr);
        auto start_atom_idx = bond->getBeginAtomIdx();
        auto end_atom_idx = bond->getEndAtomIdx();
        auto bond_item = std::make_shared<TestBondItem>(
            bond, *atom_items[start_atom_idx], *atom_items[end_atom_idx],
            test_scene->m_bond_item_settings);
        bond_item->setPos(scene_bond_items.at(idx++)->pos());
        bond_items.push_back(bond_item);
    }
    return {bond_items, test_scene};
}

std::tuple<std::shared_ptr<TestBondItem>, std::shared_ptr<TestScene>, QLineF>
createBondItem()
{
    auto test_scene = std::make_shared<TestScene>();
    test_scene->importText("CC", Format::SMILES);
    std::vector<AtomItem*> atom_items;
    for (auto* item : test_scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_items.push_back(atom_item);
        }
    }
    BOOST_TEST_REQUIRE(atom_items.size() == 2);
    RDKit::Bond* bond = *(test_scene->m_mol->bonds().begin());
    BOOST_TEST_REQUIRE(bond != nullptr);
    auto bond_item = std::make_shared<TestBondItem>(
        bond, *atom_items[0], *atom_items[1], test_scene->m_bond_item_settings);
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

BOOST_AUTO_TEST_CASE(test_arePointsOnSameSideOfLine)
{
    auto [bond_item, test_scene, trimmed_line] = createBondItem();
    QPointF line_endpoint(1, 1);
    QPointF above(0, 1);
    QPointF also_above(1, 2);
    QPointF below(1, 0);
    QPointF also_below(2, 1);
    BOOST_TEST(
        bond_item->arePointsOnSameSideOfLine(above, also_above, line_endpoint));
    BOOST_TEST(
        bond_item->arePointsOnSameSideOfLine(below, also_below, line_endpoint));
    BOOST_TEST(
        !bond_item->arePointsOnSameSideOfLine(above, below, line_endpoint));
    BOOST_TEST(!bond_item->arePointsOnSameSideOfLine(also_above, also_below,
                                                     line_endpoint));
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
    for (auto bond_item : bond_items) {
        auto bond = bond_item->m_bond;
        if (bond->getBondType() != RDKit::Bond::BondType::DOUBLE &&
            bond->getBondType() != RDKit::Bond::BondType::AROMATIC) {
            continue;
        }
        const RDKit::ROMol& molecule = bond->getOwningMol();
        RDKit::RingInfo* ring_info = molecule.getRingInfo();
        auto ring_idx = bond_item->findBestRingForBond(molecule, ring_info);
        BOOST_REQUIRE(ring_idx != 2);
    }
}

} // namespace sketcher
} // namespace schrodinger
