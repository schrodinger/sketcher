#define BOOST_TEST_MODULE test_widget_utils
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <memory>

#include <QApplication>
#include <QButtonGroup>
#include <QEvent>
#include <QMenu>
#include <QScreen>
#include <QSignalSpy>
#include <QTest>

#include "../test_common.h"
#include "schrodinger/sketcher/widget/custom_monomer_button.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

std::ostream& operator<<(std::ostream& os, const QPoint& point)
{
    return os << "QPoint(" << point.x() << ", " << point.y() << ")";
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Observe popup visibility and replay suppression for both menus and widgets.
 */
class PopupStateObserver : public QObject
{
  public:
    int m_show_count = 0;
    bool m_replay_suppressed = false;

  protected:
    bool eventFilter(QObject* watched, QEvent* event) override
    {
        if (event->type() == QEvent::Show) {
            ++m_show_count;
        } else if (event->type() == QEvent::Hide) {
            m_replay_suppressed = static_cast<QWidget*>(watched)->testAttribute(
                Qt::WA_NoMouseReplay);
        }
        return QObject::eventFilter(watched, event);
    }
};

/**
 * Clicking either kind of popup button again dismisses its popup without
 * replaying the press to the button and reopening the popup on release.
 */
BOOST_DATA_TEST_CASE(popup_dismiss_by_button,
                     boost::unit_test::data::make({false, true}),
                     custom_monomer)
{
    for (const bool initially_checked : {false, true}) {
        QWidget window;
        window.resize(400, 300);
        QButtonGroup group;
        std::unique_ptr<QToolButton> button;
        QWidget* popup = nullptr;
        if (custom_monomer) {
            button = std::make_unique<CustomMonomerButton>(&window);
            popup = button->findChild<QMenu*>();
        } else {
            auto tool_button = std::make_unique<ToolButtonWithPopup>(&window);
            popup = new QWidget(tool_button.get(), Qt::Popup);
            popup->resize(180, 80);
            tool_button->setPopupWidget(popup);
            button = std::move(tool_button);
        }
        BOOST_REQUIRE(popup != nullptr);
        group.addButton(button.get());
        button->setGeometry(40, 40, 180, 30);
        button->setCheckable(true);
        button->setChecked(initially_checked);
        window.show();
        BOOST_REQUIRE(QTest::qWaitForWindowExposed(&window));
        const auto position = button->mapTo(&window, button->rect().center());
        if (!custom_monomer && !initially_checked) {
            // The first click selects a tool; subsequent clicks open its popup.
            QTest::mouseClick(window.windowHandle(), Qt::LeftButton,
                              Qt::NoModifier, position);
            BOOST_TEST(button->isChecked());
            BOOST_TEST(!popup->isVisible());
        }
        const bool checked = button->isChecked();
        QSignalSpy clicked(button.get(), &QToolButton::clicked);
        PopupStateObserver observer;
        popup->installEventFilter(&observer);

        // Window events route through the active popup and exercise Qt's
        // platform-dependent replay handling.
        for (int opening = 1; opening <= 2; ++opening) {
            QTest::mouseClick(window.windowHandle(), Qt::LeftButton,
                              Qt::NoModifier, position);
            BOOST_REQUIRE(QTest::qWaitForWindowExposed(popup));
            BOOST_REQUIRE(popup->isVisible());
            BOOST_TEST(clicked.count() == opening);
            BOOST_TEST(observer.m_show_count == opening);

            QTest::mousePress(window.windowHandle(), Qt::LeftButton,
                              Qt::NoModifier, position);
            QApplication::processEvents();
            BOOST_TEST(!popup->isVisible());
            // Offscreen Qt cannot replay presses, so also check suppression
            // at dismissal to cover the normal test platform.
            BOOST_TEST(observer.m_replay_suppressed);
            QTest::mouseRelease(window.windowHandle(), Qt::LeftButton,
                                Qt::NoModifier, position);
            QApplication::processEvents();
            BOOST_TEST(!popup->isVisible());
            BOOST_TEST(clicked.count() == opening);
            BOOST_TEST(observer.m_show_count == opening);
            BOOST_TEST(!button->isDown());
            BOOST_TEST(button->isChecked() == checked);
        }

        // A later press elsewhere must retain normal outside-click replay.
        QTest::mouseClick(window.windowHandle(), Qt::LeftButton, Qt::NoModifier,
                          position);
        BOOST_REQUIRE(QTest::qWaitForWindowExposed(popup));
        const QPoint outside_position(350, 250);
        QTest::mousePress(window.windowHandle(), Qt::LeftButton, Qt::NoModifier,
                          outside_position);
        BOOST_TEST(!popup->isVisible());
        BOOST_TEST(!observer.m_replay_suppressed);
        QTest::mouseRelease(window.windowHandle(), Qt::LeftButton,
                            Qt::NoModifier, outside_position);
    }
}

/**
 * The position calculation that ModalDialog::moveOntoScreen() used before it
 * was extracted into clamp_position_to_bounds()
 */
static QPoint previous_modal_dialog_position(const QRect& frame,
                                             const QRect& bounds)
{
    const int x =
        std::max(std::min(frame.left(), bounds.right() - frame.width() + 1),
                 bounds.left());
    const int y =
        std::max(std::min(frame.top(), bounds.bottom() - frame.height() + 1),
                 bounds.top());
    return {x, y};
}

/**
 * Clamping moves a rectangle the minimum distance to fit, and keeps the
 * top-left corner in bounds when the rectangle is too large. This matches the
 * calculation ModalDialog used before it was shared.
 */
BOOST_AUTO_TEST_CASE(clamp_position_to_bounds_edges)
{
    const QRect bounds(100, 50, 800, 600);
    // inside: unchanged
    BOOST_TEST(clamp_position_to_bounds(QRect(200, 100, 100, 100), bounds) ==
               QPoint(200, 100));
    // past each edge
    BOOST_TEST(clamp_position_to_bounds(QRect(50, 100, 100, 100), bounds) ==
               QPoint(100, 100));
    BOOST_TEST(clamp_position_to_bounds(QRect(200, 0, 100, 100), bounds) ==
               QPoint(200, 50));
    BOOST_TEST(clamp_position_to_bounds(QRect(850, 100, 100, 100), bounds) ==
               QPoint(800, 100));
    BOOST_TEST(clamp_position_to_bounds(QRect(200, 600, 100, 100), bounds) ==
               QPoint(200, 550));
    // larger than the bounds: the top-left wins, wherever the rectangle starts
    for (const QPoint start : {QPoint(0, 0), QPoint(500, 400)}) {
        BOOST_TEST(clamp_position_to_bounds(QRect(start, QSize(1000, 700)),
                                            bounds) == QPoint(100, 50));
    }

    for (int x : {-500, 0, 100, 150, 850, 2000}) {
        for (int y : {-500, 0, 50, 100, 600, 2000}) {
            for (const QSize size :
                 {QSize(100, 100), QSize(800, 600), QSize(1000, 700)}) {
                const QRect rect(QPoint(x, y), size);
                BOOST_TEST(clamp_position_to_bounds(rect, bounds) ==
                           previous_modal_dialog_position(rect, bounds));
            }
        }
    }
}

/**
 * Popups go below the anchor if they fit, then above, then on the side with
 * more room, and are finally clamped to the bounds
 */
BOOST_AUTO_TEST_CASE(calc_popup_position_placement)
{
    const QRect bounds(0, 0, 800, 600);
    const QSize size(200, 100);
    // below
    BOOST_TEST(calc_popup_position(size, QRect(100, 100, 30, 30), bounds) ==
               QPoint(100, 130));
    // exactly fits below
    BOOST_TEST(calc_popup_position(size, QRect(100, 470, 30, 30), bounds) ==
               QPoint(100, 500));
    // above
    BOOST_TEST(calc_popup_position(size, QRect(100, 500, 30, 30), bounds) ==
               QPoint(100, 400));
    // neither fits, with more room below: shown below and clamped to the
    // bottom edge
    const QSize tall(200, 450);
    BOOST_TEST(calc_popup_position(tall, QRect(100, 200, 30, 30), bounds) ==
               QPoint(100, 150));
    // neither fits, with more room above: shown above and clamped to the top
    // edge
    BOOST_TEST(calc_popup_position(tall, QRect(100, 350, 30, 30), bounds) ==
               QPoint(100, 0));
    // clamped at the right and left edges
    BOOST_TEST(calc_popup_position(size, QRect(700, 100, 30, 30), bounds) ==
               QPoint(600, 130));
    BOOST_TEST(calc_popup_position(size, QRect(-50, 100, 30, 30), bounds) ==
               QPoint(0, 130));
    // larger than the bounds: the top-left stays in bounds
    BOOST_TEST(calc_popup_position(QSize(1000, 700), QRect(300, 300, 30, 30),
                                   bounds) == QPoint(0, 0));
    // bounds that don't start at the origin, as for a second screen
    const QRect offset_bounds(1000, 200, 800, 600);
    BOOST_TEST(calc_popup_position(size, QRect(1700, 750, 30, 30),
                                   offset_bounds) == QPoint(1600, 650));
}

/**
 * Press-and-hold popups are shown next to their button using the shared
 * placement rules and the button screen's available geometry
 */
BOOST_AUTO_TEST_CASE(tool_button_popup_position)
{
    QWidget window;
    window.resize(300, 100);
    auto* button = new ToolButtonWithPopup(&window);
    button->setGeometry(40, 40, 180, 30);
    // Frameless, so that the offscreen platform's window frame margins don't
    // offset the popup's position on its first show
    auto* popup = new QWidget(button, Qt::Popup | Qt::FramelessWindowHint);
    popup->resize(180, 80);
    button->setPopupWidget(popup);
    window.show();
    BOOST_REQUIRE(QTest::qWaitForWindowExposed(&window));
    const auto bounds = button->screen()->availableGeometry();

    for (const bool near_bottom : {false, true}) {
        // put the button near the bottom of the screen so that the popup
        // doesn't fit below it
        window.move(bounds.left(), near_bottom
                                       ? bounds.bottom() - window.height() + 1
                                       : bounds.top());
        QApplication::processEvents();
        const QRect anchor(button->mapToGlobal(QPoint(0, 0)), button->size());
        QTest::mousePress(button, Qt::LeftButton);
        BOOST_REQUIRE(
            QTest::qWaitFor([popup]() { return popup->isVisible(); }, 2000));
        QTest::mouseRelease(button, Qt::LeftButton);
        BOOST_TEST(popup->pos() ==
                   calc_popup_position(popup->size(), anchor, bounds));
        BOOST_TEST(popup->pos().y() == (near_bottom
                                            ? anchor.top() - popup->height()
                                            : anchor.bottom() + 1));
        BOOST_TEST(popup->pos().x() == anchor.left());
        popup->hide();
    }
}

} // namespace sketcher
} // namespace schrodinger
