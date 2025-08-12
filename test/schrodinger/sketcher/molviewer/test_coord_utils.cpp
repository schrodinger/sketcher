#define BOOST_TEST_MODULE Test_Sketcher

#include <cmath>

#include <rdkit/GraphMol/Depictor/RDDepictor.h>

#include <QPointF>

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
// Boost doesn't know how to print QPoints
BOOST_TEST_DONT_PRINT_LOG_VALUE(QPointF);

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_are_points_on_same_side_of_line)
{
    QPointF line_endpoint(1, 1);
    QPointF above(0, 1);
    QPointF also_above(1, 2);
    QPointF below(1, 0);
    QPointF also_below(2, 1);
    BOOST_TEST(
        are_points_on_same_side_of_line(above, also_above, line_endpoint));
    BOOST_TEST(
        are_points_on_same_side_of_line(below, also_below, line_endpoint));
    BOOST_TEST(!are_points_on_same_side_of_line(above, below, line_endpoint));
    BOOST_TEST(!are_points_on_same_side_of_line(also_above, also_below,
                                                line_endpoint));
}

BOOST_AUTO_TEST_CASE(test_rotate_conformer_radians,
                     *boost::unit_test::tolerance(0.0001))
{
    auto mol = rdkit_extensions::to_rdkit("C");
    RDDepict::compute2DCoords(*mol);
    auto& conf = mol->getConformer();
    auto& coord = conf.getAtomPos(0);
    coord = RDGeom::Point3D(0.0, 2.0, 0.0);
    RDGeom::Point3D center(0.0, 1.0, 0.0);
    rotate_conformer_radians(M_PI / 2.0, center, conf);

    BOOST_TEST(coord.x == -1.0);
    BOOST_TEST(coord.y == 1.0);
    BOOST_TEST(coord.z == 0.0);
}

/**
 * test that when points are given at very similar angles around the origin, the
 * best_placing_around_origin returns a point in the middle of the first angle
 * interval in the list, even if it's not precisely the biggest
 */
BOOST_AUTO_TEST_CASE(test_best_placing_around_origin_similar_angles,
                     *boost::unit_test::tolerance(0.001))
{
    QLineF line(QPointF(0, 0), QPointF(0, 1));
    /* test small variations approximating 120 degrees intervals. 60 degrees
     should always be returned as midpoint of the first interval*/
    std::vector<std::vector<float>> all_angles = {{0, 120, 240},
                                                  {0, 120.1, 240.1},
                                                  {0, 119.9, 240.1},
                                                  {0, 119.9, 239.9},
                                                  {0, 120.1, 239.9}};
    for (const auto& angles : all_angles) {
        std::vector<QPointF> points;
        for (auto angle : angles) {
            line.setAngle(angle);
            points.push_back(line.p2());
        }
        line.setP2(best_placing_around_origin(points));
        BOOST_TEST(line.angle() == 60);
    }
} // namespace sketcher

} // namespace sketcher
} // namespace schrodinger
