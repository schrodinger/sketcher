#define BOOST_TEST_MODULE test_widget_utils
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>

#include <QButtonGroup>
#include <QEvent>
#include <QMenu>
#include <QSignalSpy>
#include <QTest>

#include "../test_common.h"
#include "schrodinger/sketcher/widget/custom_monomer_button.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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

} // namespace sketcher
} // namespace schrodinger
