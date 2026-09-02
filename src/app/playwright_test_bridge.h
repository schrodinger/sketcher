/* -------------------------------------------------------------------------
 * Playwright-only test bridge for the WebAssembly Sketcher application.
 *
 * These functions expose Qt object-name based interaction to browser tests.
 * They are intentionally isolated from application behavior; production UI
 * code neither calls nor depends on this interface.
 * ---------------------------------------------------------------------------
 */

#pragma once

#include <string>

namespace schrodinger::sketcher
{
class SketcherWidget;

namespace playwright_test_bridge
{
void click_button(SketcherWidget& sketcher, const std::string& object_name);
void send_mouse_press(SketcherWidget& sketcher, const std::string& object_name);
void close_active_popups();
std::string get_widget_rect(SketcherWidget& sketcher,
                            const std::string& object_name);
std::string get_widget_state(SketcherWidget& sketcher,
                             const std::string& object_name);
void set_widget_text(SketcherWidget& sketcher, const std::string& object_name,
                     const std::string& text);
std::string get_menu_action_rect(SketcherWidget& sketcher,
                                 const std::string& object_name_or_text);
std::string get_atom_rect(SketcherWidget& sketcher, int atom_index);
std::string get_bond_rect(SketcherWidget& sketcher, int bond_index);
std::string get_rendered_atom_geometry(SketcherWidget& sketcher);
std::string get_rendered_bond_geometry(SketcherWidget& sketcher);
std::string visible_widget_names(SketcherWidget& sketcher);
std::string clipboard_text();
void set_clipboard_text(const std::string& text);
} // namespace playwright_test_bridge
} // namespace schrodinger::sketcher
