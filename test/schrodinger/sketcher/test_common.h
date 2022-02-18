#pragma once

#include <fstream>

#include <QtGlobal>

#include <RDGeneral/RDLog.h>

#include "schrodinger/test/testfiles.h"
#include "test_markers.h"

/* Use this fixture to construct global fixture which would provide
 * QApplication for a whole suite of a given test file.
 * It should always be used in the BOOST_GLOBAL_FIXTURE
 * This fixture is not suitable as a global fixture if test depends on
 * QApplication.
 */
class Test_Sketcher_global_fixture : public Test_QAPP_DisplayRequiredFixture
{
  public:
    Test_Sketcher_global_fixture() : Test_QAPP_DisplayRequiredFixture(true) {}
};

class SilenceRDKitLogging
{
  public:
    SilenceRDKitLogging()
    {
        if (rdErrorLog == nullptr) {
            RDLog::InitLogs();
        }
        boost::logging::disable_logs("rdApp.*");
    };
    ~SilenceRDKitLogging() { boost::logging::enable_logs("rdApp.*"); };
};

std::string read_testfile(const std::string& filename)
{
    std::ifstream fh(schrodinger::test::mmshare_testfile(filename));
    return std::string(std::istreambuf_iterator<char>(fh),
                       std::istreambuf_iterator<char>());
}
