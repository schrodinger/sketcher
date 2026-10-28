#define BOOST_TEST_MODULE test_blankable_spin_box
#include <boost/test/unit_test.hpp>

#include <memory>

#include <QProxyStyle>
#include <QStyleFactory>
#include <QStyleOptionSpinBox>
#include <QTest>
#include <QToolButton>

#include "../test_common.h"
#include "schrodinger/sketcher/widget/blankable_spin_box.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Supply deterministic stacked or side-by-side spin box button geometry.
 */
class SpinBoxStyle : public QProxyStyle
{
  public:
    SpinBoxStyle(bool side_by_side) :
        QProxyStyle(QStyleFactory::create("Fusion")),
        m_side_by_side(side_by_side)
    {
    }

    QRect subControlRect(ComplexControl control, const QStyleOptionComplex* opt,
                         SubControl sub_control,
                         const QWidget* widget = nullptr) const override
    {
        auto* spin_opt = qstyleoption_cast<const QStyleOptionSpinBox*>(opt);
        if (control != CC_SpinBox || spin_opt == nullptr) {
            return QProxyStyle::subControlRect(control, opt, sub_control,
                                               widget);
        }
        auto edit_rect = opt->rect.adjusted(2, 2, -2, -2);
        if (spin_opt->buttonSymbols == QAbstractSpinBox::NoButtons) {
            return sub_control == SC_SpinBoxEditField ? edit_rect : QRect();
        }
        auto button_width = 20;
        auto arrows_width = m_side_by_side ? 2 * button_width : button_width;
        auto arrows_left = edit_rect.right() + 1 - arrows_width;
        if (sub_control == SC_SpinBoxEditField) {
            edit_rect.setRight(arrows_left - 1);
            return edit_rect;
        }
        if (sub_control == SC_SpinBoxUp || sub_control == SC_SpinBoxDown) {
            auto down = sub_control == SC_SpinBoxDown;
            auto height =
                m_side_by_side ? edit_rect.height() : edit_rect.height() / 2;
            auto x = arrows_left + (m_side_by_side && down ? button_width : 0);
            auto y = edit_rect.top() + (!m_side_by_side && down ? height : 0);
            return QRect(x, y, button_width, height);
        }
        return QProxyStyle::subControlRect(control, opt, sub_control, widget);
    }

  private:
    bool m_side_by_side;
};

class TestSpinBox : public BlankableSpinBox
{
  public:
    using BlankableSpinBox::initStyleOption;
    using BlankableSpinBox::m_clear_btn;
};

} // namespace

/**
 * Both arrows must remain clickable, and clearing must still blank the value.
 */
static void check_clear_button(QStyle& style)
{
    for (auto symbols :
         {QAbstractSpinBox::UpDownArrows, QAbstractSpinBox::PlusMinus,
          QAbstractSpinBox::NoButtons}) {
        TestSpinBox spin_box;
        spin_box.setStyle(&style);
        spin_box.setButtonSymbols(symbols);
        spin_box.setRange(-10, 10);
        spin_box.setOptionalValue(0);
        spin_box.show();
        for (auto width : {150, 250}) {
            spin_box.resize(width, 40);
            QApplication::processEvents();
            auto* clear_btn = spin_box.m_clear_btn;
            BOOST_REQUIRE(clear_btn->isVisible());
            QStyleOptionSpinBox opt;
            spin_box.initStyleOption(&opt);
            auto edit_rect =
                style.subControlRect(QStyle::CC_SpinBox, &opt,
                                     QStyle::SC_SpinBoxEditField, &spin_box);
            BOOST_CHECK(edit_rect.contains(clear_btn->geometry()));
            if (symbols != QAbstractSpinBox::NoButtons) {
                for (auto sub_control :
                     {QStyle::SC_SpinBoxUp, QStyle::SC_SpinBoxDown}) {
                    auto arrow_rect = style.subControlRect(
                        QStyle::CC_SpinBox, &opt, sub_control, &spin_box);
                    BOOST_CHECK(!clear_btn->geometry().intersects(arrow_rect));
                    BOOST_CHECK(spin_box.childAt(arrow_rect.center()) ==
                                nullptr);
                    QTest::mouseClick(&spin_box, Qt::LeftButton, Qt::NoModifier,
                                      arrow_rect.center());
                    BOOST_CHECK(spin_box.optionalValue() ==
                                (sub_control == QStyle::SC_SpinBoxUp ? 1 : 0));
                }
            }
        }
        QTest::mouseClick(spin_box.m_clear_btn, Qt::LeftButton);
        BOOST_CHECK(!spin_box.optionalValue().has_value());
        BOOST_CHECK(!spin_box.m_clear_btn->isVisible());
    }
}

BOOST_AUTO_TEST_CASE(test_arrow_layouts)
{
    for (auto side_by_side : {false, true}) {
        BOOST_TEST_CONTEXT("side_by_side=" << side_by_side)
        {
            SpinBoxStyle style(side_by_side);
            check_clear_button(style);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_available_styles)
{
    for (const auto& name : QStyleFactory::keys()) {
        BOOST_TEST_CONTEXT("style=" << name.toStdString())
        {
            std::unique_ptr<QStyle> style(QStyleFactory::create(name));
            BOOST_REQUIRE(style != nullptr);
            check_clear_button(*style);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
