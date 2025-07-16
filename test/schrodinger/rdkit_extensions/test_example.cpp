// @copyright Schrodinger, LLC - All Rights Reserved

#define BOOST_TEST_MODULE test_example

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>

#include "schrodinger/rdkit_extensions/example.h"

using namespace schrodinger::rdkit_extensions;

BOOST_DATA_TEST_CASE(test_dependency_test,
                     boost::unit_test::data::make({"c1ccccc1", "C1=CC=CC=C1"}))
{
    BOOST_TEST(dependency_test(sample));

    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(0, 1, 0);
    Eigen::Vector3d c(1, 0, 0);
    BOOST_TEST((a - b).cross(a - c).norm() != 0.0);

    int d = 1;
    int e = 2;
    BOOST_TEST(d != e);
}
