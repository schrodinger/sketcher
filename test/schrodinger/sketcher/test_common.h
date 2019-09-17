#ifndef TEST_COMMON_H
#define TEST_COMMON_H
#include <QtGlobal>
#include "test_markers.h"

#define CREATE_QAPP                                                            \
    Test_QAPP_DisplayRequiredFixture obj(false);                               \
    if (!obj.hasDisplay()) {                                                   \
        return;                                                                \
    }

#endif // TEST_COMMON_H
