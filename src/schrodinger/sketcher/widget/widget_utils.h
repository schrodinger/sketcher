#pragma once

#include <QComboBox>

#include "schrodinger/sketcher/definitions.h"

class QAbstractButton;
class QButtonGroup;

namespace schrodinger
{
namespace sketcher
{

/**
 * Convenience function for either checking a specific button or unchecking all
 * buttons in a group.
 *
 * Specifically, supply a button pointer and a button group. If the button
 * pointer is `nullptr`, then uncheck every button in the group. Otherwise,
 * check the button.
 *
 * @param button A button belonging to `group`, or `nullptr`
 * @param group A button group containing `button`
 */
void SKETCHER_API check_button_or_uncheck_group(QAbstractButton* button,
                                                QButtonGroup* group);

/**
 * Clear a button's hover and pressed state (but not checked state) when its
 * asynchronous popup hides by adding an event filter to the popup. Without this
 * event filter, if the, e.g., "+ Custom" menu closed due to the Sketcher window
 * losing focus, the blue background would remain on the button until the next
 * time the user moused over it.
 *
 * Note that the popup owns the event filter.
 */
void SKETCHER_API clear_button_state_on_popup_hide(QAbstractButton* button,
                                                   QWidget* popup);

/**
 * Switch the combo box to the entry with the specified Qt.UserRole data
 *
 * Note that this method is intentionally defined in the header so that the
 * template can be specialized in other files.
 */
template <typename T>
static void set_combo_box_data(QComboBox* const combo, T data)
{
    auto variant = QVariant::fromValue(data);
    combo->setCurrentIndex(combo->findData(variant));
};

} // namespace sketcher
} // namespace schrodinger
