#define BOOST_TEST_MODULE Test_Sketcher

#include <QCoreApplication>
#include <QCursor>
#include <QMouseEvent>
#include <QRegularExpression>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/resizable_model_dialog.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(QRect);
BOOST_TEST_DONT_PRINT_LOG_VALUE(QSize);

namespace schrodinger
{
namespace sketcher
{

class FramelessDialog : public ResizableModelDialog
{
  public:
    FramelessDialog() :
        ResizableModelDialog(nullptr, Qt::Dialog | Qt::FramelessWindowHint)
    {
        setAttribute(Qt::WA_DeleteOnClose, false);
        layout()->setSizeConstraint(QLayout::SetNoConstraint);
        setMinimumSize(100, 80);
        setMaximumSize(500, 400);
        enableBorderResizing();
        // Contents created after the handles must not cover the resize border.
        auto* contents = new QWidget(this);
        contents->setGeometry(0, 0, 320, 240);
        contents->setCursor(Qt::CrossCursor);
        setGeometry(100, 120, 320, 240);
        show();
        QCoreApplication::processEvents();
    }
};

static QWidget* resize_handle(FramelessDialog& dialog, Qt::Edges edges)
{
    auto* handle = dialog.findChild<QWidget*>(
        "dialog_resize_handle_" + QString::number(static_cast<int>(edges)));
    BOOST_REQUIRE(handle != nullptr);
    return handle;
}

static void send_mouse_event(QWidget* widget, QEvent::Type type,
                             QPoint global_position, Qt::MouseButton button,
                             Qt::MouseButtons buttons)
{
    QMouseEvent event(type, widget->mapFromGlobal(global_position),
                      global_position, button, buttons, Qt::NoModifier);
    QCoreApplication::sendEvent(widget, &event);
}

static void drag_handle(QWidget* handle, QPoint delta)
{
    const auto start = handle->mapToGlobal(handle->rect().center());
    send_mouse_event(handle, QEvent::MouseButtonPress, start, Qt::LeftButton,
                     Qt::LeftButton);
    send_mouse_event(handle, QEvent::MouseMove, start + delta, Qt::NoButton,
                     Qt::LeftButton);
    send_mouse_event(handle, QEvent::MouseButtonRelease, start + delta,
                     Qt::LeftButton, Qt::NoButton);
}

/**
 * Every edge and corner exposes its cursor and keeps the opposite edge fixed.
 */
BOOST_AUTO_TEST_CASE(resize_edges_and_corners)
{
    FramelessDialog dialog;
    struct TestCase {
        Qt::Edges edges;
        Qt::CursorShape cursor;
        QRect expected;
    };
    const TestCase cases[] = {
        {Qt::LeftEdge, Qt::SizeHorCursor, {130, 120, 290, 240}},
        {Qt::RightEdge, Qt::SizeHorCursor, {100, 120, 350, 240}},
        {Qt::TopEdge, Qt::SizeVerCursor, {100, 140, 320, 220}},
        {Qt::BottomEdge, Qt::SizeVerCursor, {100, 120, 320, 260}},
        {Qt::LeftEdge | Qt::TopEdge, Qt::SizeFDiagCursor, {130, 140, 290, 220}},
        {Qt::RightEdge | Qt::TopEdge,
         Qt::SizeBDiagCursor,
         {100, 140, 350, 220}},
        {Qt::LeftEdge | Qt::BottomEdge,
         Qt::SizeBDiagCursor,
         {130, 120, 290, 260}},
        {Qt::RightEdge | Qt::BottomEdge,
         Qt::SizeFDiagCursor,
         {100, 120, 350, 260}}};

    for (const auto& test_case : cases) {
        dialog.setGeometry(100, 120, 320, 240);
        auto* handle = resize_handle(dialog, test_case.edges);
        const auto center = handle->mapTo(&dialog, handle->rect().center());
        BOOST_TEST(dialog.childAt(center) == handle);
        BOOST_TEST(handle->cursor().shape() == test_case.cursor);
        drag_handle(handle, {30, 20});
        BOOST_TEST(dialog.geometry() == test_case.expected);
        // A further movement after release must not continue the drag.
        send_mouse_event(handle, QEvent::MouseMove,
                         handle->mapToGlobal(QPoint(1000, 1000)), Qt::NoButton,
                         Qt::NoButton);
        BOOST_TEST(dialog.geometry() == test_case.expected);
    }
    BOOST_TEST(dialog.childAt(160, 120)->cursor().shape() == Qt::CrossCursor);
}

/**
 * Clamp large drags without moving the opposite corner, then allow reversal.
 */
BOOST_AUTO_TEST_CASE(resize_respects_size_limits)
{
    FramelessDialog dialog;
    const Qt::Edges corners[] = {
        Qt::LeftEdge | Qt::TopEdge, Qt::RightEdge | Qt::TopEdge,
        Qt::LeftEdge | Qt::BottomEdge, Qt::RightEdge | Qt::BottomEdge};
    for (const auto edges : corners) {
        const QRect original(100, 120, 320, 240);
        dialog.setGeometry(original);
        auto* handle = resize_handle(dialog, edges);
        const auto start = handle->mapToGlobal(handle->rect().center());
        const QPoint shrink(edges.testFlag(Qt::LeftEdge) ? 1000 : -1000,
                            edges.testFlag(Qt::TopEdge) ? 1000 : -1000);
        send_mouse_event(handle, QEvent::MouseButtonPress, start,
                         Qt::LeftButton, Qt::LeftButton);
        for (const auto delta : {shrink, -shrink, QPoint()}) {
            send_mouse_event(handle, QEvent::MouseMove, start + delta,
                             Qt::NoButton, Qt::LeftButton);
            const auto expected_size = delta == shrink    ? dialog.minimumSize()
                                       : delta == -shrink ? dialog.maximumSize()
                                                          : original.size();
            BOOST_TEST(dialog.size() == expected_size);
            BOOST_TEST((edges.testFlag(Qt::LeftEdge)
                            ? dialog.geometry().right() == original.right()
                            : dialog.geometry().left() == original.left()));
            BOOST_TEST((edges.testFlag(Qt::TopEdge)
                            ? dialog.geometry().bottom() == original.bottom()
                            : dialog.geometry().top() == original.top()));
        }
        send_mouse_event(handle, QEvent::MouseButtonRelease, start,
                         Qt::LeftButton, Qt::NoButton);
    }
}

BOOST_AUTO_TEST_CASE(right_button_does_not_resize)
{
    FramelessDialog dialog;
    auto* handle = resize_handle(dialog, Qt::RightEdge | Qt::BottomEdge);
    const auto original = dialog.geometry();
    const auto start = handle->mapToGlobal(handle->rect().center());
    send_mouse_event(handle, QEvent::MouseButtonPress, start, Qt::RightButton,
                     Qt::RightButton);
    send_mouse_event(handle, QEvent::MouseMove, start + QPoint(30, 20),
                     Qt::NoButton, Qt::RightButton);
    send_mouse_event(handle, QEvent::MouseButtonRelease, start + QPoint(30, 20),
                     Qt::RightButton, Qt::NoButton);
    BOOST_TEST(dialog.geometry() == original);
}

#ifndef __EMSCRIPTEN__
/**
 * Ensure that ResizableModelDialog doesn't add any resize handles in non-WASM
 * builds
 */
BOOST_AUTO_TEST_CASE(native_dialog_uses_native_frame)
{
    ResizableModelDialog dialog;
    BOOST_TEST(!dialog.windowFlags().testFlag(Qt::FramelessWindowHint));
    BOOST_TEST(dialog
                   .findChildren<QWidget*>(
                       QRegularExpression("dialog_resize_handle_.*"))
                   .isEmpty());
}
#endif

} // namespace sketcher
} // namespace schrodinger
