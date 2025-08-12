#pragma once
#include <optional>
#include <QSpinBox>
#include "schrodinger/sketcher/definitions.h"

class QToolButton;

namespace schrodinger
{
namespace sketcher
{

/**
 * A spin box where a blank value is considered an acceptable input. This spin
 * box displays a clear button when it is non-blank, just like a QLineEdit with
 * clearButtonEnabled set.
 */
class SKETCHER_API BlankableSpinBox : public QSpinBox
{
  public:
    BlankableSpinBox(QWidget* parent = nullptr);

    // reimplement QSpinBox methods to account for the blank sentinel value.
    // Note that these methods are *not* virtual in QSpinBox, so do not call
    // them on a BlankableSpinBox pointer that has been cast to a QSpinBox
    // pointer.
    void setMinimum(int minimum);
    void setRange(int minimum, int maximum);

    /**
     * Get the current spin box value
     * @return The value enterred in the spin box, or nullopt if the spin box is
     * blank
     */
    std::optional<int> optionalValue() const;

    /**
     * Set the spin box value
     * @param optional_value The value to set the spin box to, or nullopt to
     * clear the spin box
     */
    void setOptionalValue(const std::optional<int>& optional_value);

  protected:
    QToolButton* m_clear_btn = nullptr;

    // prevent clients from calling the value method, as they should use
    // optionalValue instead
    using QSpinBox::value;

    /**
     * Display the clear button if and only if the spin box contains a non-empty
     * value
     */
    void updateClearButtonVisibility();

    // override QSpinBox methods
    void resizeEvent(QResizeEvent* event) override;
    QString textFromValue(int value) const override;
    int valueFromText(const QString& text) const override;
    void stepBy(int steps) override;
    QValidator::State validate(QString& input, int& pos) const override;
    QAbstractSpinBox::StepEnabled stepEnabled() const override;
};

/**
 * A spin box that applies the same right margin as BlankableSpinBox (which
 * makes space for the clear button), but without allowing blank values or
 * displaying a clear button. UnblankableSpinBox works just like a regular
 * QSpinBox, but it aligns properly when placed above or below a
 * BlankableSpinBox.
 */
class SKETCHER_API UnblankableSpinBox : public QSpinBox
{
  public:
    UnblankableSpinBox(QWidget* parent = nullptr);
};

} // namespace sketcher
} // namespace schrodinger
