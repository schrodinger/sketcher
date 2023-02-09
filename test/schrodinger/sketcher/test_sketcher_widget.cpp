#define BOOST_TEST_MODULE sketcher_widget_test

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"

#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "test_common.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

class TestSketcherWidget : public SketcherWidget
{
  public:
    TestSketcherWidget() : SketcherWidget(){};
    using SketcherWidget::m_scene;
    using SketcherWidget::m_watermark_item;
};

/**
 * Verify that we can produce an instance of the class without error.
 */
BOOST_AUTO_TEST_CASE(test_constructor)
{
    schrodinger::sketcher::SketcherWidget sketcher_wdg;
}

/**
 * Verify that watermark visibility is trigged on inclusion/removal of atoms
 */
BOOST_AUTO_TEST_CASE(test_watermark)
{
    TestSketcherWidget sk;
    // Without the event loop, we need to manually trigger Scene::changed
    QList<QRectF> region;
    sk.m_scene->changed(region);
    // sketcher starts out empty
    BOOST_TEST(sk.m_watermark_item->isVisible());
    sk.m_scene->importText("c1ccncc1", Format::SMILES);
    sk.m_scene->changed(region);
    BOOST_TEST(!sk.m_watermark_item->isVisible());
    sk.m_scene->clearInteractiveItems();
    sk.m_scene->changed(region);
    BOOST_TEST(sk.m_watermark_item->isVisible());
}