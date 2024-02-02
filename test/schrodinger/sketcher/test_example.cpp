#define BOOST_TEST_MODULE test_example

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/example.h"

using namespace schrodinger::sketcher;

BOOST_AUTO_TEST_CASE(test_boost)
{
    BOOST_TEST(boost_library_works());
}

BOOST_AUTO_TEST_CASE(test_fmt)
{
    BOOST_TEST(fmt_library_works());
}

BOOST_AUTO_TEST_CASE(test_maeparser)
{
    // BOOST_TEST(maeparser_library_works());
}

BOOST_AUTO_TEST_CASE(test_qt6)
{
    // BOOST_TEST(qt6_library_works());
}

BOOST_DATA_TEST_CASE(test_rdkit,
                     boost::unit_test::data::make({"c1ccccc1", "C1=CC=CC=C1"}))
{
    // BOOST_TEST(rdkit_library_works());
}

BOOST_AUTO_TEST_CASE(test_zstd)
{
    BOOST_TEST(zstd_library_works());
}