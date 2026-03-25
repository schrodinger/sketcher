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
 * best_placing_around_origin returns the position closest to horizontal,
 * even if it's not the first interval in the list
 */
BOOST_AUTO_TEST_CASE(test_best_placing_around_origin_similar_angles,
                     *boost::unit_test::tolerance(0.001))
{
    QLineF line(QPointF(0, 0), QPointF(0, 1));
    /* test small variations approximating 120 degrees intervals.
     * With points at ~0°, ~120°, ~240°, the gaps have midpoints at ~60°, ~180°,
     * ~300° 180° is horizontal (distance 0), so it should always be returned
     */
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
        BOOST_TEST(line.angle() == 180);
    }
}

/**
 * test that when multiple angle intervals are tied (equal size), the one
 * closest to horizontal (0° or 180°) is preferred
 */
BOOST_AUTO_TEST_CASE(test_best_placing_around_origin_prefers_horizontal,
                     *boost::unit_test::tolerance(0.001))
{
    QLineF line(QPointF(0, 0), QPointF(0, 1));

    // Test configuration: points at 0°, 90°, 180°, 270°
    // This creates four equal 90° gaps with midpoints at 45°, 135°, 225°, 315°
    // All midpoints are equidistant (45°) from horizontal, so first one wins
    // Gap 0: 0° to 90° → midpoint 45° (dist=45) ← first, should win
    // Gap 1: 90° to 180° → midpoint 135° (dist=45)
    // Gap 2: 180° to 270° → midpoint 225° (dist=45)
    // Gap 3: 270° to 360° → midpoint 315° (dist=45)
    std::vector<float> angles = {0, 90, 180, 270};
    std::vector<QPointF> points;
    for (auto angle : angles) {
        line.setAngle(angle);
        points.push_back(line.p2());
    }

    line.setP2(best_placing_around_origin(points));
    BOOST_TEST(line.angle() == 45.0);

    // Test configuration: points at 45°, 135°, 225°, 315°
    // After sorting, angles = {45, 135, 225, 315, 405}
    // Gap 0: 45° to 135° → midpoint 90° (dist=90)
    // Gap 1: 135° to 225° → midpoint 180° (dist=0, horizontal) ← should win
    // Gap 2: 225° to 315° → midpoint 270° (dist=90)
    // Gap 3: 315° to 405° → midpoint 360°=0° (dist=0, horizontal)
    // Gaps 1 and 3 are tied, Gap 1 comes first
    angles = {45, 135, 225, 315};
    points.clear();
    for (auto angle : angles) {
        line.setAngle(angle);
        points.push_back(line.p2());
    }

    line.setP2(best_placing_around_origin(points));
    BOOST_TEST(line.angle() == 180.0);

    // Test configuration: points at 0°, 180°
    // This creates two equal 180° gaps with midpoints at 90°, 270°
    // Gap 0: 0° to 180° → midpoint 90° (dist=90) ← first, should win
    // Gap 1: 180° to 360° → midpoint 270° (dist=90)
    // Both equidistant from horizontal, so first one wins
    angles = {0, 180};
    points.clear();
    for (auto angle : angles) {
        line.setAngle(angle);
        points.push_back(line.p2());
    }

    line.setP2(best_placing_around_origin(points));
    BOOST_TEST(line.angle() == 90.0);
} // namespace sketcher

} // namespace sketcher
} // namespace schrodinger
