#define BOOST_TEST_MODULE sgroup_item

#include <QtGlobal>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/ROMol.h>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestSGroupItem : public SGroupItem
{
  public:
    using SGroupItem::m_field_data_text;
    using SGroupItem::m_label;
    using SGroupItem::m_repeat;

    TestSGroupItem(const RDKit::SubstanceGroup& sgroup, const Fonts& fonts,
                   QGraphicsItem* parent = nullptr) :
        SGroupItem(sgroup, fonts, parent)
    {
    }
};

BOOST_AUTO_TEST_CASE(test_SRU_COP_labels)
{
    auto test_scene = TestScene::getScene();
    auto smiles = "CCCCC |Sg:n:1::ht,Sg:n:2:5-8:hh,Sg:co:3::eu|";
    auto mol = rdkit_extensions::to_rdkit(smiles);
    RDDepict::compute2DCoords(*mol);
    auto sgroups = RDKit::getSubstanceGroups(*mol);
    BOOST_REQUIRE(sgroups.size() == 3);
    {
        auto sgroup_item = TestSGroupItem(sgroups[0], test_scene->m_fonts);
        BOOST_TEST(sgroup_item.m_label.toStdString() == "n"); // ht
        BOOST_TEST(sgroup_item.m_repeat.toStdString() == "");
        BOOST_TEST(sgroup_item.m_field_data_text.toStdString() == "");
    }
    {
        auto sgroup_item = TestSGroupItem(sgroups[1], test_scene->m_fonts);
        BOOST_TEST(sgroup_item.m_label.toStdString() == "5-8");
        BOOST_TEST(sgroup_item.m_repeat.toStdString() == "hh");
        BOOST_TEST(sgroup_item.m_field_data_text.toStdString() == "");
    }
    {
        auto sgroup_item = TestSGroupItem(sgroups[2], test_scene->m_fonts);
        BOOST_TEST(sgroup_item.m_label.toStdString() == "co");
        BOOST_TEST(sgroup_item.m_repeat.toStdString() == "eu");
        BOOST_TEST(sgroup_item.m_field_data_text.toStdString() == "");
    }
}

BOOST_AUTO_TEST_CASE(test_DAT_labels)
{
    auto test_scene = TestScene::getScene();
    std::string molblock{R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -24.4583 19.7917 0 0
M  V30 2 C -23.1247 20.5617 0 0
M  V30 3 C -21.791 19.7917 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 3) FIELDNAME=value FIELDDISP="  -21.7910   19.7917    -
M  V30 DAU   ALL  0       0" FIELDDATA="foo"
M  V30 END SGROUP
M  V30 END CTAB
M  END
$$$$
)MDL"};
    auto mol = rdkit_extensions::to_rdkit(molblock);
    auto sgroups = RDKit::getSubstanceGroups(*mol);
    BOOST_REQUIRE(sgroups.size() == 1);
    auto sgroup_item = TestSGroupItem(sgroups[0], test_scene->m_fonts);
    BOOST_TEST(sgroup_item.m_label.toStdString() == "");
    BOOST_TEST(sgroup_item.m_repeat.toStdString() == "");
    BOOST_TEST(sgroup_item.m_field_data_text.toStdString() == "foo");
}

} // namespace sketcher
} // namespace schrodinger
