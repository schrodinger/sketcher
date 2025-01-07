#define BOOST_TEST_MODULE Test_Sketcher

#include <string>
#include <utility>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <QMarginsF>
#include <QRectF>
#include <QString>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
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
    using AtomItem::getQueryLabel;
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

    TestAtomItem(RDKit::Atom* atom, Fonts& fonts,
                 const AtomDisplaySettings& settings,
                 QGraphicsItem* parent = nullptr) :
        AtomItem(atom, fonts, settings, parent)
    {
    }
};

std::pair<std::vector<std::shared_ptr<TestAtomItem>>,
          std::shared_ptr<TestScene>>
createAtomItems(std::string smiles)
{
    auto test_scene = TestScene::getScene();
    import_mol_text(test_scene->m_mol_model, smiles);
    std::vector<std::shared_ptr<TestAtomItem>> atom_items;
    for (auto atom : test_scene->m_mol_model->getMol()->atoms()) {
        BOOST_TEST_REQUIRE(atom != nullptr);

        auto atom_item = std::make_shared<TestAtomItem>(
            atom, test_scene->m_fonts,
            *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
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
    BOOST_TEST(atom_item->m_charge_and_radical_label_text == "(2•) 2+");
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

BOOST_AUTO_TEST_CASE(test_updateCachedData_atomLabel)
{
    for (const auto input_str : {"* |$Aryl_p$|", "[#0] |$Aryl_p$|"}) {
        auto [atom_item, test_scene] = createAtomItem(input_str);
        BOOST_TEST(atom_item->m_main_label_text == "Aryl_p");
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
}

BOOST_AUTO_TEST_CASE(test_findPositionInEmptySpace,
                     *boost::unit_test::tolerance(0.01))
{
    auto scene = TestScene::getScene();
    import_mol_text(scene->m_mol_model, "CCC", Format::SMILES);
    std::vector<AtomItem*> atom_items;
    for (auto* item : scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_items.push_back(atom_item);
        }
    }

    auto central_atom = atom_items.at(1);
    auto first_atom_position = atom_items.at(0)->pos();
    auto second_atom_position = atom_items.at(2)->pos();
    auto empty_space_position = central_atom->findPositionInEmptySpace(false);

    for (const auto& pos : {first_atom_position, second_atom_position}) {
        BOOST_TEST(QLineF(pos, empty_space_position).length() == 71.9);
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

BOOST_AUTO_TEST_CASE(test_enhanced_chirality_labels)
{
    // test that enhanced chirality labels are correctly set, replacing the
    // chirality label for OR and AND but not for ABS
    std::string smiles =
        "C[C@H](N)[C@H](C)[C@H](C)[C@@H](N)[C@H](C)N |a:1,3,o3:7,9,&1:5|";
    auto [atom_items, scene] = createAtomItems(smiles);
    std::vector<std::string> labels = {"", "(S)",  "", "(R)",  "", "and 1",
                                       "", "or 3", "", "or 3", "", ""};

    for (unsigned int i = 0; i < atom_items.size(); ++i) {
        BOOST_TEST(atom_items.at(i)->m_chirality_label_text.toStdString() ==
                   labels.at(i));
    }
}

BOOST_AUTO_TEST_CASE(test_get_query_label)
{
    /*check that setting an advanced property on an allowed list switches the
     * query label to a smarts*/
    auto props = std::make_shared<AtomQueryProperties>();
    props->query_type = QueryType::ALLOWED_LIST;
    props->allowed_list = {Element::C, Element::N};

    std::vector<std::string> labels = {"[C,N]", "\"[#6,#7;X2]\""};
    for (unsigned int i = 0; i < labels.size(); ++i) {
        if (i == 1) {
            props->num_connections = 2;
        }
        auto test_scene = TestScene::getScene();
        auto [atom, _enh_stereo] = create_atom_with_properties(props);
        std::shared_ptr<RDKit::Atom::QUERYATOM_QUERY> query;
        query.reset(atom->getQuery()->copy());
        test_scene->m_mol_model->addAtom(query, {0, 0, 0});
        auto model_atom = test_scene->m_mol_model->getMol()->getAtomWithIdx(0);
        auto atom_item = std::make_shared<TestAtomItem>(
            const_cast<RDKit::Atom*>(model_atom), test_scene->m_fonts,
            *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
        BOOST_TEST(atom_item->getQueryLabel().toStdString() == labels.at(i));
    }
}

} // namespace sketcher
} // namespace schrodinger
