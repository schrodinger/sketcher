#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/data/test_case.hpp>

#include <QUndoStack>
#include <QShowEvent>
#include <QWheelEvent>
#include <QtWidgets/qtestsupport_widgets.h>

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/view.h"
// #include "schrodinger/sketcher/molviewer/scene.h"
// #include "schrodinger/sketcher/model/mol_model.h"
// #include "schrodinger/sketcher/model/sketcher_model.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestView : public View
{
  public:
    TestView(QGraphicsScene* scene, QWidget* parent = nullptr) :
        View(scene, parent)
    {
    }
    using View::m_delayed_fit_to_screen;
    using View::m_initial_geometry_set;
    using View::showEvent;
    using View::wheelEvent;
};

/**
 * When calling fitToScreen() before the View has been shown, make sure that a
 * delayed fitToScreen() call is scheduled for immediately after the View is
 * shown.
 */
BOOST_AUTO_TEST_CASE(test_delayed_fit_to_screen)
{
    auto scene = TestScene::getScene();
    TestView view(scene.get());
    BOOST_TEST(!view.m_delayed_fit_to_screen);
    BOOST_TEST(!view.m_initial_geometry_set);
    view.fitToScreen();
    BOOST_TEST(view.m_delayed_fit_to_screen);
    BOOST_TEST(!view.m_initial_geometry_set);

    QShowEvent show_event;
    view.showEvent(&show_event);
    BOOST_TEST(!view.m_delayed_fit_to_screen);
    BOOST_TEST(view.m_initial_geometry_set);
}

/**
 * SKETCH-2772: a wheel zoom must keep the scene point under the cursor fixed,
 * so users can focus on a region by zooming while pointing at it.
 */
BOOST_AUTO_TEST_CASE(test_wheel_zoom_keeps_cursor_point_fixed)
{
    auto scene = TestScene::getScene();
    TestView view(scene.get());
    view.resize(400, 400);
    view.show();
    [[maybe_unused]] auto exposed = QTest::qWaitForWindowExposed(&view);
    // The View defaults to an empty sceneRect (real usage sets one via
    // fitToScreen). Cursor-anchored zoom uses sceneRect to position the
    // viewport, so prime a non-empty one for the test.
    view.setSceneRect(QRectF(-200, -200, 400, 400));

    auto cursor = QPointF(300, 100); // offset from the view center (200, 200)
    auto scene_pos_before = view.mapToScene(cursor.toPoint());

    // Negative deltaY = wheel down = zoom out, which avoids the zoom-in cap.
    auto global_cursor = view.viewport()->mapToGlobal(cursor.toPoint());
    QWheelEvent wheel_event(cursor, global_cursor, QPoint(0, -120),
                            QPoint(0, -120), Qt::NoButton, Qt::NoModifier,
                            Qt::NoScrollPhase, false);
    view.wheelEvent(&wheel_event);

    auto scene_pos_after = view.mapToScene(cursor.toPoint());
    // Qt rounds scroll positions to integers, so the correction is accurate
    // to within ~1 pixel in scene coordinates (~1% relative error).
    BOOST_CHECK_CLOSE(scene_pos_after.x(), scene_pos_before.x(), 1.0);
    BOOST_CHECK_CLOSE(scene_pos_after.y(), scene_pos_before.y(), 1.0);
}

} // namespace sketcher
} // namespace schrodinger
