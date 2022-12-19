#define BOOST_TEST_MODULE sketcher_widget_test
#include <boost/test/unit_test.hpp>
#include "test_common.h"
#include "schrodinger/sketcher/sketcher_widget.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

/**
 * Verify that we can produce an instance of the class without error.
 */
BOOST_AUTO_TEST_CASE(constructor)
{

    schrodinger::sketcher::SketcherWidget sketcher_wdg;
}
