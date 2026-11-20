#define BOOST_TEST_MODULE test_multi_select_combo_box
#include <boost/test/unit_test.hpp>

#include <memory>
#include <utility>
#include <vector>

#include <QAbstractItemView>
#include <QApplication>
#include <QFontMetrics>
#include <QImage>
#include <QMouseEvent>
#include <QScrollBar>
#include <QSignalBlocker>
#include <QSignalSpy>
#include <QStandardItemModel>
#include <QStyle>
#include <QStyleFactory>
#include <QStyleOptionComboBox>
#include <QTest>

#include "../test_common.h"
#include "schrodinger/sketcher/widget/multi_select_combo_box.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * A visible multi-select combo box with three items, whose data is 1, 2 and
 * 4
 */
struct ComboFixture {
    ComboFixture()
    {
        combo.addItem("Peptide", 1);
        combo.addItem("Nucleic acid", 2);
        combo.addItem("Structure", 4);
        combo.show();
    }

    ~ComboFixture()
    {
        combo.hidePopup();
    }

    /**
     * Click the given item in the open dropdown
     */
    void clickItem(int row)
    {
        auto* view = combo.view();
        BOOST_REQUIRE(view->isVisible());
        QTest::mouseClick(
            view->viewport(), Qt::LeftButton, Qt::NoModifier,
            view->visualRect(combo.model()->index(row, 0)).center());
        QApplication::processEvents();
    }

    /**
     * Move the mouse over the given item in the open dropdown
     */
    void hoverItem(int row)
    {
        auto* view = combo.view();
        BOOST_REQUIRE(view->isVisible());
        auto* viewport = view->viewport();
        const auto pos =
            view->visualRect(combo.model()->index(row, 0)).center();
        QMouseEvent move(QEvent::MouseMove, pos, viewport->mapToGlobal(pos),
                         Qt::NoButton, Qt::NoButton, Qt::NoModifier);
        QApplication::sendEvent(viewport, &move);
    }

    /**
     * Press the given key in the open dropdown
     */
    void pressKey(Qt::Key key, Qt::KeyboardModifiers modifiers = {})
    {
        auto* view = combo.view();
        BOOST_REQUIRE(view->isVisible());
        QTest::keyClick(view, key, modifiers);
    }

    /**
     * @return The rows that the dropdown highlights
     */
    std::vector<int> getHighlightedRows()
    {
        std::vector<int> rows;
        for (const auto& index :
             combo.view()->selectionModel()->selectedIndexes()) {
            rows.push_back(index.row());
        }
        return rows;
    }

    /**
     * @return The dropdown's current row, which keyboard input applies to
     */
    int getCurrentRow()
    {
        return combo.view()->currentIndex().row();
    }

    std::vector<bool> getCheckStates()
    {
        std::vector<bool> states;
        for (int row = 0; row < combo.count(); ++row) {
            states.push_back(combo.isItemChecked(row));
        }
        return states;
    }

    MultiSelectComboBox combo;
};

/**
 * Added items are checkable and unchecked, and the display text joins the
 * checked items in item order, or reads "Any" when all are checked or "None"
 * when none are
 */
BOOST_AUTO_TEST_CASE(display_text)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    for (int row = 0; row < combo.count(); ++row) {
        auto flags = combo.model()->index(row, 0).flags();
        BOOST_TEST(flags.testFlag(Qt::ItemIsUserCheckable));
        BOOST_TEST(!combo.isItemChecked(row));
    }
    BOOST_TEST(combo.getDisplayText() == "None");
    combo.setItemChecked(2, true);
    combo.setItemChecked(0, true);
    BOOST_TEST(combo.getDisplayText() == "Peptide, Structure");
    combo.setItemChecked(1, true);
    BOOST_TEST(combo.getDisplayText() == "Any");
    combo.setItemChecked(0, false);
    BOOST_TEST(combo.getDisplayText() == "Nucleic acid, Structure");
    // the display text never comes from the current item
    BOOST_TEST(combo.currentText() == "Peptide");

    // there's room for the widest text that can be displayed
    const auto metrics = combo.fontMetrics();
    BOOST_TEST(combo.sizeHint().width() >
               metrics.horizontalAdvance("Nucleic acid, Structure"));
}

/**
 * The data API checks exactly the given items and reports them in item order,
 * emitting the change signal once per call
 */
BOOST_AUTO_TEST_CASE(checked_data)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    QSignalSpy spy(&combo, &MultiSelectComboBox::checkedItemsChanged);
    combo.setCheckedData({4, 1});
    BOOST_TEST(spy.count() == 1);
    BOOST_CHECK(combo.checkedData() == (QList<QVariant>{1, 4}));
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, false, true}),
               boost::test_tools::per_element());

    // no change, no signal
    combo.setCheckedData({1, 4});
    BOOST_TEST(spy.count() == 1);

    combo.setCheckedData({2});
    BOOST_TEST(spy.count() == 2);
    BOOST_CHECK(combo.checkedData() == (QList<QVariant>{2}));

    // programmatic single-item changes emit the signal too, unless blocked
    combo.setItemChecked(0, true);
    BOOST_TEST(spy.count() == 3);
    {
        QSignalBlocker blocker(&combo);
        combo.setItemChecked(0, false);
    }
    BOOST_TEST(spy.count() == 3);
    BOOST_TEST(!combo.isItemChecked(0));
}

/**
 * Clicking an item in the dropdown toggles it and keeps the dropdown open. A
 * press and release on different rows doesn't toggle anything, nor does a
 * click on a disabled item.
 */
BOOST_AUTO_TEST_CASE(click_toggles)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    QSignalSpy spy(&combo, &MultiSelectComboBox::checkedItemsChanged);
    combo.showPopup();
    auto* view = combo.view();
    BOOST_REQUIRE(view->isVisible());

    fixture.clickItem(1);
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{false, true, false}),
               boost::test_tools::per_element());
    BOOST_TEST(spy.count() == 1);
    BOOST_TEST(view->isVisible());
    fixture.clickItem(0);
    fixture.clickItem(1);
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, false, false}),
               boost::test_tools::per_element());
    BOOST_TEST(spy.count() == 3);
    BOOST_TEST(view->isVisible());

    // press on one row and release on another
    auto center = [&](int row) {
        return view->visualRect(combo.model()->index(row, 0)).center();
    };
    QTest::mousePress(view->viewport(), Qt::LeftButton, Qt::NoModifier,
                      center(1));
    QTest::mouseRelease(view->viewport(), Qt::LeftButton, Qt::NoModifier,
                        center(2));
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, false, false}),
               boost::test_tools::per_element());
    BOOST_TEST(view->isVisible());

    // a disabled item can't be toggled
    auto* model = qobject_cast<QStandardItemModel*>(combo.model());
    BOOST_REQUIRE(model != nullptr);
    model->item(2)->setEnabled(false);
    fixture.clickItem(2);
    BOOST_TEST(!combo.isItemChecked(2));
    BOOST_TEST(view->isVisible());
}

/**
 * Space toggles the current item, keeps the dropdown open and highlights the
 * current item, so that the next Up or Down moves the highlight. Escape
 * closes the dropdown.
 */
BOOST_AUTO_TEST_CASE(space_toggles)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.showPopup();
    auto* view = combo.view();
    BOOST_REQUIRE(view->isVisible());
    BOOST_TEST(fixture.getHighlightedRows().empty());

    fixture.pressKey(Qt::Key_Space);
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, false, false}),
               boost::test_tools::per_element());
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{0},
               boost::test_tools::per_element());
    BOOST_TEST(view->isVisible());

    fixture.pressKey(Qt::Key_Down);
    BOOST_TEST(fixture.getCurrentRow() == 1);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{1},
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Space);
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, true, false}),
               boost::test_tools::per_element());
    BOOST_TEST(view->isVisible());

    fixture.pressKey(Qt::Key_Escape);
    BOOST_TEST(!view->isVisible());
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{true, true, false}),
               boost::test_tools::per_element());
}

/**
 * The dropdown opens with no highlight. The first Up or Down highlights the
 * current item without moving it, and later ones move the highlight.
 */
BOOST_AUTO_TEST_CASE(arrows_move_highlight)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.showPopup();
    BOOST_TEST(fixture.getHighlightedRows().empty());
    BOOST_TEST(fixture.getCurrentRow() == 0);

    fixture.pressKey(Qt::Key_Down);
    BOOST_TEST(fixture.getCurrentRow() == 0);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{0},
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Down);
    fixture.pressKey(Qt::Key_Down);
    BOOST_TEST(fixture.getCurrentRow() == 2);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{2},
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Up);
    BOOST_TEST(fixture.getCurrentRow() == 1);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{1},
               boost::test_tools::per_element());
    // navigating doesn't check anything
    BOOST_TEST(combo.checkedData().isEmpty());

    // reopening keeps the current item but hides the highlight until the next
    // key press, and the keypad's arrows behave like the main ones
    combo.hidePopup();
    combo.showPopup();
    BOOST_TEST(fixture.getHighlightedRows().empty());
    BOOST_TEST(fixture.getCurrentRow() == 1);
    fixture.pressKey(Qt::Key_Up, Qt::KeypadModifier);
    BOOST_TEST(fixture.getCurrentRow() == 1);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{1},
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Up, Qt::KeypadModifier);
    BOOST_TEST(fixture.getCurrentRow() == 0);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{0},
               boost::test_tools::per_element());
}

/**
 * Hovering over an item highlights it and makes it current, so that Space
 * toggles it and the first Up or Down moves away from it
 */
BOOST_AUTO_TEST_CASE(hover_highlights)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.showPopup();
    fixture.hoverItem(2);
    BOOST_TEST(fixture.getCurrentRow() == 2);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{2},
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Space);
    BOOST_TEST(fixture.getCheckStates() ==
                   (std::vector<bool>{false, false, true}),
               boost::test_tools::per_element());
    fixture.pressKey(Qt::Key_Up);
    BOOST_TEST(fixture.getCurrentRow() == 1);
    BOOST_TEST(fixture.getHighlightedRows() == std::vector<int>{1},
               boost::test_tools::per_element());

    // a hover highlight from before the dropdown closed doesn't reappear
    combo.hidePopup();
    combo.showPopup();
    BOOST_TEST(fixture.getHighlightedRows().empty());
}

/**
 * Enter and Return toggle the current item and close the dropdown, unless the
 * item is disabled
 */
BOOST_AUTO_TEST_CASE(enter_toggles_and_closes)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    auto* view = combo.view();
    combo.showPopup();
    fixture.pressKey(Qt::Key_Down);
    fixture.pressKey(Qt::Key_Down);
    BOOST_REQUIRE(fixture.getCurrentRow() == 1);
    std::vector<bool> expected{false, false, false};
    for (const auto& [key, modifiers] :
         std::vector<std::pair<Qt::Key, Qt::KeyboardModifiers>>{
             {Qt::Key_Enter, Qt::NoModifier},
             {Qt::Key_Return, Qt::NoModifier},
             {Qt::Key_Enter, Qt::KeypadModifier}}) {
        BOOST_TEST_CONTEXT("key=" << key << " modifiers=" << modifiers.toInt())
        {
            if (!view->isVisible()) {
                combo.showPopup();
            }
            fixture.pressKey(key, modifiers);
            expected[1] = !expected[1];
            BOOST_TEST(fixture.getCheckStates() == expected,
                       boost::test_tools::per_element());
            BOOST_TEST(!view->isVisible());
        }
    }

    // a disabled item can't be toggled, so the dropdown stays open
    auto* model = qobject_cast<QStandardItemModel*>(combo.model());
    BOOST_REQUIRE(model != nullptr);
    model->item(1)->setEnabled(false);
    combo.showPopup();
    BOOST_REQUIRE(fixture.getCurrentRow() == 1);
    fixture.pressKey(Qt::Key_Return);
    BOOST_TEST(fixture.getCheckStates() == expected,
               boost::test_tools::per_element());
    BOOST_TEST(view->isVisible());
}

/**
 * The user can uncheck the last checked item, which leaves no items checked,
 * and can then check another item
 */
BOOST_AUTO_TEST_CASE(uncheck_last_item)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.setCheckedData({2});
    combo.showPopup();
    fixture.clickItem(1);
    BOOST_TEST(combo.checkedData().isEmpty());
    BOOST_TEST(combo.getDisplayText() == "None");
    BOOST_TEST(combo.view()->isVisible());
    fixture.clickItem(0);
    BOOST_CHECK(combo.checkedData() == (QList<QVariant>{1}));
    BOOST_TEST(combo.getDisplayText() == "Peptide");
    combo.view()->setCurrentIndex(combo.model()->index(0, 0));
    QTest::keyClick(combo.view(), Qt::Key_Space);
    BOOST_TEST(combo.checkedData().isEmpty());
    combo.hidePopup();
}

/**
 * The dropdown leaves space above its first item and below its last item, and
 * is tall enough to show every item without scrolling
 */
BOOST_AUTO_TEST_CASE(dropdown_padding)
{
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.showPopup();
    auto* view = combo.view();
    BOOST_REQUIRE(view->isVisible());
    const auto first_rect = view->visualRect(combo.model()->index(0, 0));
    const auto last_rect =
        view->visualRect(combo.model()->index(combo.count() - 1, 0));
    BOOST_TEST(first_rect.top() > 0);
    BOOST_TEST(last_rect.bottom() < view->viewport()->height() - 1);
    BOOST_TEST(view->verticalScrollBar()->maximum() == 0);
}

/**
 * With a style that draws the dropdown like a menu, such as Fusion (which the
 * WASM build uses), the space around the items has the same color as the
 * items' background
 */
BOOST_AUTO_TEST_CASE(dropdown_spacing_matches_item_background)
{
    // Declared before the fixture so that the style outlives the combo box
    std::unique_ptr<QStyle> fusion_style(QStyleFactory::create("Fusion"));
    BOOST_REQUIRE(fusion_style != nullptr);
    ComboFixture fixture;
    auto& combo = fixture.combo;
    combo.setStyle(fusion_style.get());
    QStyleOptionComboBox option;
    option.initFrom(&combo);
    BOOST_REQUIRE(
        fusion_style->styleHint(QStyle::SH_ComboBox_Popup, &option, &combo));
    combo.showPopup();
    auto* view = combo.view();
    BOOST_REQUIRE(view->isVisible());
    const auto first_rect = view->visualRect(combo.model()->index(0, 0));
    const auto second_rect = view->visualRect(combo.model()->index(1, 0));
    const auto last_rect =
        view->visualRect(combo.model()->index(combo.count() - 1, 0));
    BOOST_REQUIRE(first_rect.top() > 0);
    BOOST_REQUIRE(first_rect.left() > 0);
    BOOST_REQUIRE(second_rect.top() > first_rect.bottom() + 1);
    BOOST_REQUIRE(last_rect.bottom() < view->viewport()->height() - 1);

    // QImage::pixel() takes device pixels rather than logical pixels
    const auto image = view->viewport()->grab().toImage();
    const qreal ratio = view->viewport()->devicePixelRatioF();
    auto color_at = [&image, ratio](int x, int y) {
        return image.pixel(qRound(x * ratio), qRound(y * ratio));
    };
    // Sample the middle of each row, which is past the checkboxes. The top row
    // of an item's pixels is above its checkbox and text.
    const int viewport_width = view->viewport()->width();
    const int x = viewport_width / 2;
    const auto item_color = color_at(x, first_rect.top());
    BOOST_TEST(color_at(x, first_rect.top() - 1) == item_color);
    BOOST_TEST(color_at(x, (first_rect.bottom() + second_rect.top()) / 2) ==
               item_color);
    BOOST_TEST(color_at(x, last_rect.bottom() + 1) == item_color);
    BOOST_TEST(color_at(first_rect.left() - 1, first_rect.top()) == item_color);
    BOOST_TEST(color_at(viewport_width - 1, first_rect.top()) == item_color);
}

} // namespace sketcher
} // namespace schrodinger
