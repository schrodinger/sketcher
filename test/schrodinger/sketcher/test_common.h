#pragma once

#include <QtGlobal>

/* Use this fixture to construct global fixture which would provide
 * QApplication for a whole suite of a given test file.
 * It should always be used in the BOOST_GLOBAL_FIXTURE
 * This fixture is not suitable as a global fixture if test depends on
 * QApplication.
 */
class Test_Sketcher_global_fixture : public Test_QAPP_DisplayRequiredFixture
{
  public:
    Test_Sketcher_global_fixture() : Test_QAPP_DisplayRequiredFixture(true)
    {
    }
};
