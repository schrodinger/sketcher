#include "schrodinger/sketcher/widget/blankable_spin_box.h"

#include <fmt/format.h>

#include <QLineEdit>
#include <QString>
#include <QStyle>
#include <QStyleOptionSpinBox>
#include <QToolButton>

#include <schrodinger/sketcher/molviewer/constants.h>

namespace schrodinger
{
namespace sketcher
{

/**
 * Create and return a style sheet that adds a right margin large enough to fit
 * a clear button
 * @param style The style of the widget that this style sheet is for. This is
 * needed to fetch the size of the clear button.
 */
static QString
get_style_sheet_to_make_space_for_clear_btn(const QStyle* const style)
{
    auto space_for_icon =
        style->pixelMetric(QStyle::PM_LineEditIconSize) +
        style->pixelMetric(QStyle::PM_LayoutHorizontalSpacing);
    auto style_sheet =
        fmt::format("QLineEdit {{ padding-right: {}px; }}", space_for_icon);
    return QString::fromStdString(style_sheet);
}

BlankableSpinBox::BlankableSpinBox(QWidget* parent) : QSpinBox(parent)
{

    m_clear_btn = new QToolButton(this);
#ifdef EMSCRIPTEN
    // Adding an icon to this button in a WASM build causes the QApplication to
    // crash when the button is shown, so we use the letter X instead. See
    // SKETCH-2310 for more info.
    m_clear_btn->setText("X");
#else
    m_clear_btn->setIcon(QIcon(LINE_EDIT_CLEAR_ICON_PATH));
#endif
    m_clear_btn->setCursor(Qt::ArrowCursor);

    updateClearButtonVisibility();
    connect(m_clear_btn, &QToolButton::clicked, this, &QSpinBox::clear);
    connect(this, &QSpinBox::textChanged, this,
            &BlankableSpinBox::updateClearButtonVisibility);
    auto line_edit_style_sheet =
        get_style_sheet_to_make_space_for_clear_btn(style());
    lineEdit()->setStyleSheet(line_edit_style_sheet);
}

void BlankableSpinBox::setMinimum(int minimum)
{
    QSpinBox::setMinimum(minimum - 1);
    updateClearButtonVisibility();
}

void BlankableSpinBox::setRange(int minimum, int maximum)
{
    QSpinBox::setRange(minimum - 1, maximum);
    updateClearButtonVisibility();
}

std::optional<int> BlankableSpinBox::optionalValue() const
{
    auto value = QSpinBox::value();
    if (value == minimum()) {
        return std::nullopt;
    }
    return value;
}

void BlankableSpinBox::setOptionalValue(
    const std::optional<int>& optional_value)
{
    if (optional_value.has_value()) {
        setValue(*optional_value);
    } else {
        setValue(minimum());
    }
    updateClearButtonVisibility();
}

void BlankableSpinBox::updateClearButtonVisibility()
{
    m_clear_btn->setVisible(value() != minimum());
}

void BlankableSpinBox::resizeEvent(QResizeEvent* event)
{
    // The clear button needs to be moved to the left if the arrows are present.
    auto space_for_arrows = 0;
    if (buttonSymbols() == QAbstractSpinBox::UpDownArrows) {
        QStyleOptionSpinBox opt;
        space_for_arrows =
            style()
                ->subControlRect(QStyle::CC_SpinBox, &opt, QStyle::SC_SpinBoxUp)
                .width();
    }
    QSpinBox::resizeEvent(event);
    auto size_hint = m_clear_btn->sizeHint();
    auto frame_width = style()->pixelMetric(QStyle::PM_DefaultFrameWidth);
    auto sb_rect = rect();
    auto x =
        sb_rect.right() - frame_width - size_hint.width() - space_for_arrows;
    auto y = (sb_rect.bottom() + 1 - size_hint.height()) / 2;
    m_clear_btn->move(x, y);
}

QString BlankableSpinBox::textFromValue(int value) const
{
    if (value == minimum()) {
        return "";
    }
    return QSpinBox::textFromValue(value);
}

int BlankableSpinBox::valueFromText(const QString& text) const
{
    if (text.isEmpty()) {
        return minimum();
    }
    return QSpinBox::valueFromText(text);
}

void BlankableSpinBox::stepBy(int steps)
{
    auto min = minimum();
    if (value() == min) {
        int start_val = 0;
        if (min > 0 || maximum() < 0) {
            start_val = min;
        }
        setValue(start_val);
        if (steps >= -1 && steps <= 1) {
            // if the user wanted to do a single step, then leave the value at
            // zero
            return;
        }
    }
    QSpinBox::stepBy(steps);
}

QValidator::State BlankableSpinBox::validate(QString& input, int& pos) const
{
    if (input.isEmpty()) {
        return QValidator::Acceptable;
    }
    return QSpinBox::validate(input, pos);
}

QAbstractSpinBox::StepEnabled BlankableSpinBox::stepEnabled() const
{
    auto value = QSpinBox::value();
    auto min = minimum();
    if (value == min) {
        return QSpinBox::StepUpEnabled | QSpinBox::StepDownEnabled;
    } else if (value == min + 1) {
        return QSpinBox::StepUpEnabled;
    }
    return QSpinBox::stepEnabled();
}

UnblankableSpinBox::UnblankableSpinBox(QWidget* parent) : QSpinBox(parent)
{
    auto line_edit_style_sheet =
        get_style_sheet_to_make_space_for_clear_btn(style());
    lineEdit()->setStyleSheet(line_edit_style_sheet);
}

} // namespace sketcher
} // namespace schrodinger
