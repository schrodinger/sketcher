#pragma once

#include <QComboBox>
#include <QPoint>
#include <QRect>
#include <QSize>

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
 * Add an event filter to a button's asynchronous popup to suppress replay of
 * dismissing mouse presses over the button and clear its hover and pressed
 * state (but not checked state) when the popup hides. This prevents the popup
 * from reopening on release and handles dismissal paths that do not deliver a
 * mouse leave event to the button, such as losing window focus.
 *
 * Note that the popup owns the event filter.
 */
void SKETCHER_API configure_button_popup(QAbstractButton* button,
                                         QWidget* popup);

/**
 * Calculate the position that moves a rectangle the minimum distance needed to
 * fit inside the given bounds. If the rectangle is larger than the bounds, its
 * top-left corner is kept inside the bounds and the unavoidable overflow
 * extends past the right and bottom edges.
 *
 * @param rect The rectangle to move, e.g. a window's frame geometry
 * @param bounds The area that the rectangle should fit inside
 * @return The new top-left position for the rectangle
 */
SKETCHER_API QPoint clamp_position_to_bounds(const QRect& rect,
                                             const QRect& bounds);

/**
 * Calculate where to show a popup next to an anchor. The popup's left edge is
 * aligned with the anchor's left edge, and the popup is placed:
 *   1. below the anchor, if it fits;
 *   2. otherwise above the anchor, if it fits;
 *   3. otherwise on whichever side has more room.
 * Finally, the position is clamped to the bounds (see
 * clamp_position_to_bounds()).
 *
 * @param popup_size The size of the popup
 * @param anchor_rect The anchor's geometry, in the same coordinates as bounds
 * @param bounds The area that the popup should fit inside
 * @return The top-left position for the popup
 */
SKETCHER_API QPoint calc_popup_position(const QSize& popup_size,
                                        const QRect& anchor_rect,
                                        const QRect& bounds);

/**
 * Show a popup next to an anchor widget, positioned using
 * calc_popup_position() within the available geometry of the anchor's screen.
 * (On WASM builds, each canvas is its own screen.) This only positions the
 * popup; it never resizes it.
 *
 * @param popup The popup to show, which must be a top-level widget (e.g. a
 * Qt::Popup)
 * @param anchor The widget to show the popup next to, e.g. the button that
 * opened it
 */
SKETCHER_API void show_popup_next_to_widget(QWidget* popup, QWidget* anchor);

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
