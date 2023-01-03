#define BOOST_TEST_MODULE Test_Sketcher

#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <string>
#include <utility>

#include <QRectF>
#include <QMarginsF>
#include <QString>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);
// Boost doesn't know how to print QStrings
BOOST_TEST_DONT_PRINT_LOG_VALUE(QString);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

class TestScene : public Scene
{
  public:
    using Scene::m_atom_item_settings;
    using Scene::m_bond_item_settings;
    using Scene::m_fonts;
    using Scene::m_mol;
};

class TestAtomItem : public AtomItem
{
  public:
    using AtomItem::m_bounding_rect;
    using AtomItem::m_label_is_visible;
    using AtomItem::m_main_label_rect;
    using AtomItem::m_main_label_text;
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

std::pair<std::shared_ptr<TestAtomItem>, std::shared_ptr<TestScene>>
createAtomItem(std::string smiles)
{
    auto test_scene = std::make_shared<TestScene>();
    test_scene->importText(smiles, Format::SMILES);
    RDKit::Atom* atom = *(test_scene->m_mol->atoms().begin());
    BOOST_TEST_REQUIRE(atom != nullptr);
    auto atom_item = std::make_shared<TestAtomItem>(
        atom, test_scene->m_fonts, test_scene->m_atom_item_settings);
    return std::make_pair(atom_item, test_scene);
}

BOOST_AUTO_TEST_CASE(test_updateCachedData_LabeledN)
{
    auto [atom_item, test_scene] = createAtomItem("N");
    BOOST_TEST(atom_item->m_main_label_text == "N");
    BOOST_TEST(!atom_item->m_main_label_rect.isNull());
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
    BOOST_TEST(atom_item->m_subrects.empty());
    BOOST_TEST(!atom_item->m_label_is_visible);
    BOOST_TEST(!(atom_item->m_shape.isEmpty()));
    BOOST_TEST(!(atom_item->m_bounding_rect.isNull()));
}

} // namespace sketcher
} // namespace schrodinger
