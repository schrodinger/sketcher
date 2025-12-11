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

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
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
    using AtomItem::getSubrects;
    using AtomItem::getTooltip;
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
    using AtomItem::m_query_label_text;
    using AtomItem::m_selection_highlighting_path;
    using AtomItem::m_shape;

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
    std::map<const RDKit::Atom*, QPointF> atom_to_position;
    // get atom positions
    for (auto* item : test_scene->items()) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_to_position[atom_item->getAtom()] = atom_item->pos();
        }
    }

    std::vector<std::shared_ptr<TestAtomItem>> atom_items;
    for (auto atom : test_scene->m_mol_model->getMol()->atoms()) {
        BOOST_TEST_REQUIRE(atom != nullptr);

        auto atom_item = std::make_shared<TestAtomItem>(
            atom, test_scene->m_fonts,
            *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
        atom_item->setPos(atom_to_position[atom]);
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

std::pair<std::shared_ptr<TestAtomItem>, std::shared_ptr<TestScene>>
createAtomItem(std::shared_ptr<AtomQueryProperties> props)
{
    auto test_scene = TestScene::getScene();
    auto [atom, _enh_stereo] = create_atom_with_properties(props);
    std::shared_ptr<RDKit::Atom::QUERYATOM_QUERY> query;
    query.reset(atom->getQuery()->copy());
    test_scene->m_mol_model->addAtom(query, {0, 0, 0});
    auto model_atom = test_scene->m_mol_model->getMol()->getAtomWithIdx(0);
    auto atom_item = std::make_shared<TestAtomItem>(
        const_cast<RDKit::Atom*>(model_atom), test_scene->m_fonts,
        *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
    return std::make_pair(atom_item, test_scene);
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

    BOOST_TEST(!atom_item->getSubrects().empty());
    BOOST_TEST(atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));

    // QPainterPath.contains(QRect) will return false if the rect touches the
    // edge of the path, so we have to shrink all of our rectangles by one pixel
    // to ensure that they're completely inside the path
    QMarginsF margins(1, 1, 1, 1);
    auto& shape = atom_item->m_shape;
    BOOST_TEST(shape.contains(atom_item->m_main_label_rect - margins));
    for (QRectF subrect : atom_item->getSubrects()) {
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
    BOOST_TEST(atom_item->getSubrects().empty());
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
    BOOST_TEST(!atom_item->getSubrects().empty());
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

    BOOST_TEST(!atom_item->getSubrects().empty());
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
    BOOST_TEST(!atom_item->getSubrects().empty());
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
        BOOST_TEST(!atom_item->getSubrects().empty());
        BOOST_TEST(atom_item->m_label_is_visible);
        BOOST_TEST(!(atom_item->m_shape.isEmpty()));
        BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));
    }
}

BOOST_AUTO_TEST_CASE(test_deuterium_and_tritium_display)
{

    /*check the deuterium and tritium are rendered with D and T when the
     * relevant flag is on*/

    auto [atom_items, test_scene] = createAtomItems("[2H]-[3H]");

    BOOST_TEST(atom_items.at(0)->m_main_label_text == "H");
    BOOST_TEST(atom_items.at(0)->m_isotope_label_text == "2");

    BOOST_TEST(atom_items.at(1)->m_main_label_text == "H");
    BOOST_TEST(atom_items.at(1)->m_isotope_label_text == "3");

    auto display_settings(
        *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
    display_settings.m_show_symbol_for_H_isotopes = true;
    test_scene->m_sketcher_model->setAtomDisplaySettings(display_settings);

    atom_items[0]->updateCachedData();
    atom_items[1]->updateCachedData();

    BOOST_TEST(atom_items[0]->m_main_label_text == "D");
    BOOST_TEST(atom_items.at(0)->m_isotope_label_text.toStdString() == "");
    BOOST_TEST(atom_items.at(0)->m_isotope_label_rect.isNull());

    BOOST_TEST(atom_items[1]->m_main_label_text == "T");
    BOOST_TEST(atom_items.at(1)->m_isotope_label_text.toStdString() == "");
    BOOST_TEST(atom_items.at(1)->m_isotope_label_rect.isNull());
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
        auto [atom_item, test_scene] = createAtomItem(props);

        auto label = atom_item->getQueryLabel().toStdString();
        // on some platforms the label could be elided, let's check up to the
        // 6th character
        BOOST_TEST(label.substr(0, 6) == labels.at(i).substr(0, 6));
    }
}

BOOST_AUTO_TEST_CASE(test_no_Hs_on_queries)
{
    /*check that setting a query results in an atomic label that shows the
     * element without any implicit Hs*/
    auto [atom_item, test_scene] = createAtomItem("[#6+]");
    BOOST_TEST(atom_item->m_main_label_text == "C");
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
    BOOST_TEST(!atom_item->m_charge_and_radical_label_rect.isNull());
    BOOST_TEST(atom_item->m_charge_and_radical_label_text == "+");
    BOOST_TEST(atom_item->m_H_count_label_rect.isNull());
    BOOST_TEST(atom_item->m_H_label_rect.isNull());
}

BOOST_AUTO_TEST_CASE(test_queries_rendering)
{
    /*check the rendering of different types of queries*/

    { // specific element query with advanced properties
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::SPECIFIC_ELEMENT;
        props->element = Element::C;
        props->ring_count_type = QueryCount::EXACTLY;
        props->ring_count_exact_val = 2;
        props->num_connections = 2;
        props->charge = 1;
        auto [atom_item, test_scene] = createAtomItem(props);
        BOOST_TEST(atom_item->m_main_label_text == "C");
        BOOST_TEST(atom_item->m_query_label_text == "X2&R2");
        BOOST_TEST(atom_item->m_charge_and_radical_label_text == "+");
    }

    { // atom list
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::ALLOWED_LIST;
        props->allowed_list = {Element::C, Element::N};
        auto [atom_item, test_scene] = createAtomItem(props);
        BOOST_TEST(atom_item->m_main_label_text == "*");
        BOOST_TEST(atom_item->m_query_label_text == "[C,N]");
    }

    { // atom list with charge
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::ALLOWED_LIST;
        props->allowed_list = {Element::C, Element::N};
        props->charge = 1;
        auto [atom_item, test_scene] = createAtomItem(props);
        BOOST_TEST(atom_item->m_main_label_text == "*");
        BOOST_TEST(atom_item->m_query_label_text == "\"[#6,#7;+]\"");
    }

    { // atom list with enhanced stereo (SKETCH-2487)
        // Stereo properties shouldn't prevent displaying element symbols
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::ALLOWED_LIST;
        props->allowed_list = {Element::N, Element::C, Element::O, Element::F};
        props->enhanced_stereo =
            EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0);
        auto [atom_item, test_scene] = createAtomItem(props);
        BOOST_TEST(atom_item->m_main_label_text == "*");
        // Should use element symbols, not atomic numbers
        BOOST_TEST(atom_item->m_query_label_text == "[N,C,O,F]");
    }

    { // wildcard query with advanced properties
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::WILDCARD;
        props->wildcard = AtomQuery::XH;
        props->ring_count_type = QueryCount::EXACTLY;
        props->ring_count_exact_val = 4;
        props->num_connections = 2;
        props->charge = 3;
        auto [atom_item, test_scene] = createAtomItem(props);
        BOOST_TEST(atom_item->m_main_label_text == "XH");
        BOOST_TEST(atom_item->m_query_label_text == "X2&R4");
        BOOST_TEST(atom_item->m_charge_and_radical_label_text == "3+");
    }
}

BOOST_AUTO_TEST_CASE(test_chirality_label_positioning)
{
    /*check the chirality labels are positioned correctly inside the ring
    SKETCH-2576*/

    std::string smiles = "FC1(S)C(F)(S)C(F)(S)C(F)(S)C(F)(S)C1(F)S";

    auto [atom_items, scene] = createAtomItems(smiles);
    // get the positions from the scene atoms

    // get all the atoms with a chirality label (i.e. the C atoms)
    std::vector<std::shared_ptr<TestAtomItem>> carbon_items;
    for (auto item : atom_items) {
        if (item->m_chirality_label_rect.isValid()) {
            carbon_items.push_back(item);
        }
    }
    BOOST_TEST(carbon_items.size() == 6);
    // find the centroid of the C atoms
    QPointF centroid(0, 0);
    for (auto item : carbon_items) {
        centroid += item->pos();
    }
    centroid /= carbon_items.size();
    // for each C atom, check that the center of the chirality label is closer
    // to the centroid than the C atom position, i.e. the labels points towards
    // the center of the ring
    for (auto item : carbon_items) {
        auto chirality_label_center = item->m_chirality_label_rect.center();
        BOOST_TEST(QLineF(centroid, chirality_label_center).length() <
                   QLineF(centroid, item->pos()).length());
    }
}

BOOST_AUTO_TEST_CASE(test_chirality_label_tooltip)
{
    // Test that atoms with chirality labels have tooltips
    auto [atom_items, scene] = createAtomItems("CC[C@H](C)N");

    // Find the chiral center (should be at index 2)
    std::shared_ptr<TestAtomItem> chiral_atom;
    for (auto item : atom_items) {
        if (!item->m_chirality_label_text.isEmpty()) {
            chiral_atom = item;
            break;
        }
    }

    BOOST_TEST_REQUIRE(chiral_atom != nullptr);
    BOOST_TEST(chiral_atom->m_chirality_label_text == "(S)");
    BOOST_TEST(chiral_atom->m_chirality_label_rect.isValid());

    // Tooltip should show the chirality label with Stereo prefix
    BOOST_TEST(chiral_atom->toolTip() == "Stereo: (S)");
}

BOOST_AUTO_TEST_CASE(test_wildcard_atom_tooltips)
{
    // Test that wildcard atoms have descriptive tooltips
    struct TestCase {
        AtomQuery wildcard;
        QString expected_tooltip;
    };

    std::vector<TestCase> test_cases = {
        {AtomQuery::Q, "Query: Any heteroatom"},
        {AtomQuery::A, "Query: Any heavy atom"},
        {AtomQuery::X, "Query: Any halogen"},
        {AtomQuery::M, "Query: Any metal"},
        {AtomQuery::QH, "Query: Any heteroatom or hydrogen"},
        {AtomQuery::AH, "Query: Any heavy atom or hydrogen"},
        {AtomQuery::XH, "Query: Any halogen or hydrogen"},
        {AtomQuery::MH, "Query: Any metal or hydrogen"}};

    for (const auto& [wildcard, expected_tooltip] : test_cases) {
        auto props = std::make_shared<AtomQueryProperties>();
        props->query_type = QueryType::WILDCARD;
        props->wildcard = wildcard;

        auto [atom_item, test_scene] = createAtomItem(props);

        // Wildcards should have descriptive tooltips
        BOOST_TEST(atom_item->toolTip().toStdString() ==
                   expected_tooltip.toStdString());
    }
}

BOOST_AUTO_TEST_CASE(test_query_atom_with_constraints_tooltip)
{
    // Test that query atoms with additional constraints show the constraint
    // text
    auto props = std::make_shared<AtomQueryProperties>();
    props->query_type = QueryType::WILDCARD;
    props->wildcard = AtomQuery::X;
    props->ring_count_type = QueryCount::EXACTLY;
    props->ring_count_exact_val = 2;
    props->num_connections = 2;

    auto [atom_item, test_scene] = createAtomItem(props);

    BOOST_TEST(atom_item->m_main_label_text == "X");
    BOOST_TEST(atom_item->m_query_label_text == "X2&R2");

    // Should show the query constraint text with Query prefix
    BOOST_TEST(atom_item->toolTip() == "Query: X2&R2");
}

BOOST_AUTO_TEST_CASE(test_atom_list_tooltip)
{
    // Test that atom lists show the list as tooltip
    auto props = std::make_shared<AtomQueryProperties>();
    props->query_type = QueryType::ALLOWED_LIST;
    props->allowed_list = {Element::C, Element::N};

    auto [atom_item, test_scene] = createAtomItem(props);

    BOOST_TEST(atom_item->m_main_label_text == "*");
    BOOST_TEST(atom_item->m_query_label_text == "[C,N]");

    // Should show the atom list with Query prefix
    BOOST_TEST(atom_item->toolTip() == "Query: [C,N]");
}

BOOST_AUTO_TEST_CASE(test_no_tooltip_for_normal_atoms)
{
    // Test that normal atom labels don't have tooltips
    std::vector<std::string> normal_atoms = {"C", "N", "O", "S", "P"};

    for (const auto& smiles : normal_atoms) {
        auto [atom_item, scene] = createAtomItem(smiles);

        // Normal atoms should not have tooltips
        BOOST_TEST(atom_item->toolTip().isEmpty());
    }
}

BOOST_AUTO_TEST_CASE(test_chirality_takes_precedence_over_query)
{
    // Test that atoms with only chirality (no query) show only stereo tooltip
    auto [atom_items, scene] = createAtomItems("C[C@H](N)S");

    for (auto item : atom_items) {
        if (!item->m_chirality_label_text.isEmpty()) {
            // Should show chirality with Stereo prefix
            BOOST_TEST(item->toolTip() ==
                       "Stereo: " + item->m_chirality_label_text);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_atom_with_both_chirality_and_query)
{
    // Test that atoms with both chirality and query features show both in
    // tooltip
    auto props = std::make_shared<AtomQueryProperties>();
    props->query_type = QueryType::SPECIFIC_ELEMENT;
    props->element = Element::C;
    props->num_connections = 2;

    auto [atom_item, test_scene] = createAtomItem(props);

    // Manually set chirality label to simulate a chiral query atom
    atom_item->m_chirality_label_text = "(R)";
    atom_item->m_query_label_text = "X2";

    // Tooltip should show both stereo and query on separate lines
    QString tooltip = atom_item->getTooltip();
    BOOST_TEST(tooltip.contains("Stereo: (R)"));
    BOOST_TEST(tooltip.contains("Query: X2"));
}

BOOST_AUTO_TEST_CASE(test_atom_index_display)
{
    auto [atom_items, test_scene] = createAtomItems("CCN");

    // By default, atom indices should not be shown (carbons are unlabeled)
    BOOST_TEST(!atom_items[0]->m_label_is_visible); // C
    BOOST_TEST(!atom_items[1]->m_label_is_visible); // C
    BOOST_TEST(atom_items[2]->m_label_is_visible);  // N
    BOOST_TEST(atom_items[2]->m_main_label_text == "N");

    // Enable atom index display
    auto display_settings(
        *test_scene->m_sketcher_model->getAtomDisplaySettingsPtr());
    display_settings.m_show_atom_indices = true;
    test_scene->m_sketcher_model->setAtomDisplaySettings(display_settings);

    // Update all atom items
    for (auto atom_item : atom_items) {
        atom_item->updateCachedData();
    }

    // Verify that all atoms now show labels with indices appended
    BOOST_TEST(atom_items[0]->m_label_is_visible);
    BOOST_TEST(atom_items[0]->m_main_label_text == "C:0");

    BOOST_TEST(atom_items[1]->m_label_is_visible);
    BOOST_TEST(atom_items[1]->m_main_label_text == "C:1");

    BOOST_TEST(atom_items[2]->m_label_is_visible);
    BOOST_TEST(atom_items[2]->m_main_label_text == "N:2");
}

} // namespace sketcher
} // namespace schrodinger
