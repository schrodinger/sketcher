#define BOOST_TEST_MODULE Test_Sketcher

#include <string>
#include <utility>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <QMarginsF>
#include <QRectF>
#include <QString>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);
// Boost doesn't know how to print QStrings
BOOST_TEST_DONT_PRINT_LOG_VALUE(QString);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::HsDirection);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

class TestAtomItem : public AtomItem
{
  public:
    using AtomItem::determineValenceErrorIsVisible;
    using AtomItem::findHsDirection;
    using AtomItem::m_bounding_rect;
    using AtomItem::m_charge_and_radical_label_rect;
    using AtomItem::m_charge_and_radical_label_text;
    using AtomItem::m_chirality_label_rect;
    using AtomItem::m_chirality_label_text;
    using AtomItem::m_H_count_label_rect;
    using AtomItem::m_H_count_label_text;
    using AtomItem::m_H_label_rect;
    using AtomItem::m_isotope_label_rect;
    using AtomItem::m_isotope_label_text;
    using AtomItem::m_label_is_visible;
    using AtomItem::m_main_label_rect;
    using AtomItem::m_main_label_text;
    using AtomItem::m_mapping_label_rect;
    using AtomItem::m_mapping_label_text;
    using AtomItem::m_predictive_highlighting_path;
    using AtomItem::m_selection_highlighting_path;
    using AtomItem::m_shape;
    using AtomItem::m_subrects;

    TestAtomItem(RDKit::Atom* atom, Fonts& fonts, AtomItemSettings& settings,
                 QGraphicsItem* parent = nullptr) :
        AtomItem(atom, fonts, settings, parent)
    {
    }
};

std::pair<std::vector<std::shared_ptr<TestAtomItem>>,
          std::shared_ptr<TestScene>>
createAtomItems(std::string smiles)
{
    auto test_scene = std::make_shared<TestScene>();
    test_scene->m_mol_model->importFromText(smiles, Format::SMILES);
    std::vector<std::shared_ptr<TestAtomItem>> atom_items;
    for (auto atom : test_scene->m_mol_model->getMol()->atoms()) {
        BOOST_TEST_REQUIRE(atom != nullptr);
        auto atom_item = std::make_shared<TestAtomItem>(
            atom, test_scene->m_fonts, test_scene->m_atom_item_settings);
        atom_items.push_back(atom_item);
    }
    return std::make_pair(atom_items, test_scene);
}

std::pair<std::shared_ptr<TestAtomItem>, std::shared_ptr<TestScene>>
createAtomItem(std::string smiles)
{
    auto [atom_items, test_scene] = createAtomItems(smiles);
    return {atom_items.front(), test_scene};
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_LabeledN)
{
    auto [atom_item, test_scene] = createAtomItem("N");
    BOOST_TEST(atom_item->m_main_label_text == "N");
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
    BOOST_TEST(atom_item->m_isotope_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(!atom_item->m_H_count_label_rect.isNull());
    BOOST_TEST(!atom_item->m_H_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_count_label_text == "3");

    BOOST_TEST(!atom_item->m_subrects.empty());
    BOOST_TEST(atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));

    // QPainterPath.contains(QRect) will return false if the rect touches the
    // edge of the path, so we have to shrink all of our rectangles by one pixel
    // to ensure that they're completely inside the path
    QMarginsF margins(1, 1, 1, 1);
    auto& shape = atom_item->m_shape;
    BOOST_TEST(shape.contains(atom_item->m_main_label_rect - margins));
    for (QRectF subrect : atom_item->m_subrects) {
        BOOST_TEST(shape.contains(subrect - margins));
    }
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_UnlabeledC)
{
    auto [atom_item, test_scene] = createAtomItem("CC");
    BOOST_TEST(atom_item->m_main_label_text.isEmpty());
    BOOST_TEST(atom_item->m_main_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_count_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_count_label_text.isEmpty());
    BOOST_TEST(atom_item->m_subrects.empty());
    BOOST_TEST(!atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_ChargeAndIsotope)
{
    auto [atom_item, test_scene] = createAtomItem("[13C++]");
    BOOST_TEST(atom_item->m_label_is_visible);
    BOOST_TEST(!atom_item->m_main_label_text.isEmpty());
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
    BOOST_TEST(!atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(!atom_item->m_isotope_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_text == "2+");
    BOOST_TEST(atom_item->m_isotope_label_text == "13");
    BOOST_TEST(!atom_item->m_subrects.empty());
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_RGroup)
{
    auto [atom_item, test_scene] = createAtomItem("* |$_R5$|");
    BOOST_TEST(atom_item->m_main_label_text == "R5");
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
    BOOST_TEST(atom_item->m_isotope_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_count_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_label_rect.isNull());

    BOOST_TEST(!atom_item->m_subrects.empty());
    BOOST_TEST(atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_atomMapping)
{
    auto [atom_item, test_scene] = createAtomItem("[CH4:1]");
    // test that the atom label is set correctly when a mapping number is
    // present, with a C and Hs and no charge or isotopes and that labels and
    // rects are correctly created
    BOOST_TEST(atom_item->m_main_label_text == "C");
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
    BOOST_TEST(atom_item->m_isotope_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(!atom_item->m_H_count_label_rect.isNull());
    BOOST_TEST(!atom_item->m_H_label_rect.isNull());
    BOOST_TEST(!atom_item->m_subrects.empty());
    BOOST_TEST(atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));

    // test the mapping portion of the label
    BOOST_TEST(!atom_item->m_mapping_label_rect.isNull());
    BOOST_TEST(atom_item->m_mapping_label_text == "[1]");
}

BOOST_AUTO_TEST_CASE(test_findPositionInEmptySpace,
                     *boost::unit_test::tolerance(0.01))
{
    TestScene scene;
    scene.m_mol_model->importFromText("CCC", Format::SMILES);
    std::vector<AtomItem*> atom_items;
    for (auto* item : scene.items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_items.push_back(atom_item);
        }
    }

    auto central_atom = atom_items.at(1);
    auto first_atom_position = atom_items.at(0)->pos();
    auto second_atom_position = atom_items.at(2)->pos();
    auto empty_space_position = central_atom->findPositionInEmptySpace(false);

    for (const auto& pos : {first_atom_position, second_atom_position}) {
        BOOST_TEST(QLineF(pos, empty_space_position).length() == 59.5);
    }
}

BOOST_AUTO_TEST_CASE(test_findHsDirection, *boost::unit_test::tolerance(0.01))
{
    auto [atom_item, scene] = createAtomItem("N");
    BOOST_TEST(atom_item->findHsDirection() == HsDirection::RIGHT);

    auto [oxygen_atom_item, scene2] = createAtomItem("O");
    BOOST_TEST(oxygen_atom_item->findHsDirection() == HsDirection::LEFT);

    auto [atom_items, scene3] = createAtomItems("NNNN");
    BOOST_TEST(atom_items.at(0)->findHsDirection() == HsDirection::LEFT);
    BOOST_TEST(atom_items.at(1)->findHsDirection() == HsDirection::UP);
    BOOST_TEST(atom_items.at(2)->findHsDirection() == HsDirection::DOWN);
    BOOST_TEST(atom_items.at(3)->findHsDirection() == HsDirection::RIGHT);
}

BOOST_AUTO_TEST_CASE(test_chirality_label)
{
    // test that the chirality are correctly set for each atom of two
    // enantiomers (non chiral atoms should have an empty label and null rect)
    std::map<std::string, std::string> chiralities = {{"C[C@H](N)S", "(R)"},
                                                      {"C[C@@H](N)S", "(S)"}};
    for (const auto& [smiles, expected_chirality] : chiralities) {
        auto [atom_items, scene] = createAtomItems(smiles);
        auto chiral_center = atom_items.at(1);
        for (auto atom_item : atom_items) {
            auto chirality =
                (atom_item == chiral_center ? expected_chirality : "");
            BOOST_TEST(atom_item->m_chirality_label_text.toStdString() ==
                       chirality);
            BOOST_TEST(atom_item->m_chirality_label_rect.isNull() ==
                       (atom_item != chiral_center));
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
