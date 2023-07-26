#define BOOST_TEST_MODULE sketcher_widget_test

#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/sketcher_widget.h"
#include "test_common.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

using namespace schrodinger::sketcher;

class TestSketcherWidget : public SketcherWidget
{
  public:
    TestSketcherWidget() : SketcherWidget(){};
};

BOOST_AUTO_TEST_CASE(test_widget)
{
    TestSketcherWidget sk;
}
