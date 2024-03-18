// @copyright Schrodinger, LLC - All Rights Reserved

#define BOOST_TEST_MODULE example_widget_test

#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/example_widget.h"
#include "qapplication_required_fixture.h"

using namespace schrodinger::sketcher;

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

class TestExampleWidget : public ExampleWidget
{
  public:
    TestExampleWidget() : ExampleWidget(){};
};

BOOST_AUTO_TEST_CASE(test_widget)
{
    TestExampleWidget widget;
}
