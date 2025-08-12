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
