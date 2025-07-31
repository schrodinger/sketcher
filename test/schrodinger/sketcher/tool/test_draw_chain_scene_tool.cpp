#define BOOST_TEST_MODULE Test_Sketcher

#include <QPointF>
#include <QList>

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/tool/draw_chain_scene_tool.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Run get_bond_chain_atom_coords and compare the output to the expected
 * coordinates.  Note that all parameters will be scaled up by
 * BOND_LENGTH * VIEW_SCALE.  In other words, a distance of 1.0 here always
 * represents the length of one bond.
 * @param start The start position for the chain
 * @param end The end position for the chain
 * @param exp_coords The expected function output
 */
void test_coords(QPointF start, QPointF end, const QList<QPointF>& exp_coords)
{
    start *= BOND_LENGTH * VIEW_SCALE;
    end *= BOND_LENGTH * VIEW_SCALE;
    auto coords = get_bond_chain_atom_coords(start, end);
    BOOST_TEST_REQUIRE(coords.length() == exp_coords.length());
    for (unsigned int i = 0; i < exp_coords.length(); ++i) {
        QPointF scaled_exp = exp_coords[i] * BOND_LENGTH * VIEW_SCALE;
        BOOST_TEST(QLineF(coords[i], scaled_exp).length() < 0.001);
    }
}

BOOST_AUTO_TEST_CASE(test_get_bond_chain_atom_coords)
{
    const qreal COS_PI_OVER_6 = qCos(M_PI / 6);
    test_coords({0, 0}, {1, 0}, {{0, 0}, {COS_PI_OVER_6, -0.5}});
    test_coords({0, 0}, {2, 0},
                {{0, 0}, {COS_PI_OVER_6, -0.5}, {2 * COS_PI_OVER_6, 0}});
    test_coords({0, 0}, {4, 0},
                {{0, 0},
                 {COS_PI_OVER_6, -0.5},
                 {2 * COS_PI_OVER_6, 0},
                 {3 * COS_PI_OVER_6, -0.5},
                 {4 * COS_PI_OVER_6, 0},
                 {5 * COS_PI_OVER_6, -0.5}});

    test_coords({0, 0}, {-2, 0},
                {{0, 0}, {-COS_PI_OVER_6, 0.5}, {-2 * COS_PI_OVER_6, 0}});
}

} // namespace sketcher
} // namespace schrodinger
