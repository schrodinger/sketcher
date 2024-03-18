// @copyright Schrodinger, LLC - All Rights Reserved

#define BOOST_TEST_MODULE test_example

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/example.h"

using namespace schrodinger::rdkit_extensions;

BOOST_DATA_TEST_CASE(test_dependency_test,
                     boost::unit_test::data::make({"c1ccccc1", "C1=CC=CC=C1"}))
{
    BOOST_TEST(dependency_test(sample));
}
