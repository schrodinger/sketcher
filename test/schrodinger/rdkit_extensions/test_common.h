#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <fstream>

#include <QtGlobal>

#include <RDGeneral/RDLog.h>

#include "schrodinger/test/testfiles.h"
#include "test_markers.h"

#define CREATE_QAPP                              \
    Test_QAPP_DisplayRequiredFixture obj(false); \
    if (!obj.hasDisplay()) {                     \
        return;                                  \
    }

#endif // TEST_COMMON_H

std::string read_testfile(const std::string& filename)
{
    std::ifstream fh(schrodinger::test::mmshare_testfile(filename));
    return std::string(std::istreambuf_iterator<char>(fh),
                       std::istreambuf_iterator<char>());
}
