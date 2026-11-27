#define BOOST_TEST_MODULE Test_Sketcher

#include <QApplication>
#include <QMouseEvent>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidget>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/custom_title_bar.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(QPoint);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::TitleBarType);

namespace schrodinger
{
namespace sketcher
{

/**
 * Drag the title bar with the left mouse button from its center by the given
 * offset
 */
static void drag_title_bar(CustomTitleBar& title_bar, const QPoint& offset)
{
    const QPointF start = title_bar.rect().center();
    const QPointF end = start + offset;
    QMouseEvent press(QEvent::MouseButtonPress, start,
                      title_bar.mapToGlobal(start), Qt::LeftButton,
                      Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&title_bar, &press);
    QMouseEvent move(QEvent::MouseMove, end, title_bar.mapToGlobal(end),
                     Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&title_bar, &move);
    QMouseEvent release(QEvent::MouseButtonRelease, end,
                        title_bar.mapToGlobal(end), Qt::LeftButton,
                        Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(&title_bar, &release);
}

/**
 * A dialog's title bar has no close button
 */
BOOST_AUTO_TEST_CASE(test_dialog_close_button)
{
    CustomTitleBar title_bar("Title", nullptr, TitleBarType::DIALOG);
    BOOST_TEST(title_bar.getCloseButton() == nullptr);
}

/**
 * A popup's title bar has a close button that fits in the title bar
 */
BOOST_AUTO_TEST_CASE(test_popup_close_button)
{
    QWidget parent;
    auto* title_bar = new CustomTitleBar("Title", &parent, TitleBarType::POPUP);
    auto* close_button = title_bar->getCloseButton();
    BOOST_REQUIRE(close_button != nullptr);
    BOOST_TEST(close_button->parent() == title_bar);
    // a multiplication sign
    BOOST_TEST(close_button->text() == QString(QChar(0x00D7)));
    // the button fits in the title bar's fixed height
    title_bar->resize(200, title_bar->height());
    parent.show();
    QApplication::processEvents();
    BOOST_TEST(title_bar->rect().contains(close_button->geometry()));
}

/**
 * A popup's title bar is shorter than a dialog's
 */
BOOST_AUTO_TEST_CASE(test_popup_height)
{
    CustomTitleBar dialog_title_bar("Title", nullptr, TitleBarType::DIALOG);
    CustomTitleBar popup_title_bar("Title", nullptr, TitleBarType::POPUP);
    BOOST_TEST(popup_title_bar.height() < dialog_title_bar.height());
}

/**
 * Dragging a dialog's title bar moves its parent, but dragging a popup's title
 * bar doesn't
 */
BOOST_DATA_TEST_CASE(test_draggable,
                     boost::unit_test::data::make({TitleBarType::DIALOG,
                                                   TitleBarType::POPUP}) ^
                         boost::unit_test::data::make({true, false}),
                     type, expect_move)
{
    QWidget parent;
    auto* layout = new QVBoxLayout(&parent);
    layout->setContentsMargins(0, 0, 0, 0);
    auto* title_bar = new CustomTitleBar("Title", &parent, type);
    layout->addWidget(title_bar);
    parent.resize(200, 100);
    parent.move(100, 100);
    parent.show();
    QApplication::processEvents();
    const QPoint offset(40, 30);

    const auto start_pos = parent.pos();
    drag_title_bar(*title_bar, offset);
    const auto expected_pos = expect_move ? start_pos + offset : start_pos;
    BOOST_TEST(parent.pos() == expected_pos);
}

} // namespace sketcher
} // namespace schrodinger
