#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_file_import_export

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

BOOST_AUTO_TEST_CASE(test_get_import_name_filters)
{
    auto filters = get_import_name_filters();
    BOOST_TEST(filters.contains(";;"));
    BOOST_TEST(filters.split(";;").size() == 9);
    BOOST_TEST(filters.contains("*.sdf"));
    BOOST_TEST(filters.contains("*.mdl"));
}
