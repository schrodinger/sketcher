#define BOOST_TEST_MODULE Test_Sketcher

#include <boost/test/data/test_case.hpp>

#include <QUndoStack>
#include <QShowEvent>

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

} // namespace sketcher
} // namespace schrodinger
