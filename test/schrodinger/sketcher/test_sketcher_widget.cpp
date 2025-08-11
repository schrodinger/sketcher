// @copyright Schrodinger, LLC - All Rights Reserved

#define BOOST_TEST_MODULE sketcher_widget_test

#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/sketcher_widget.h"
#include "qapplication_required_fixture.h"

using namespace schrodinger::sketcher;

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

class TestSketcherWidget : public SketcherWidget
{
  public:
    TestSketcherWidget() : SketcherWidget(){};
};

BOOST_AUTO_TEST_CASE(test_widget)
{
    TestSketcherWidget widget;
}
