#pragma once

namespace schrodinger
{
namespace sketcher
{

enum class SelectMode {
    SELECT,      // Shift + click
    DESELECT,    //
    TOGGLE,      // Ctrl + click
    SELECT_ONLY, // normal click (no Shift or Ctrl)
};

/**
 * Possible selection tools to equip
 */
enum class SelectionTool {
    RECTANGLE,
    ELLIPSE,
    LASSO,
    FRAGMENT,
};

} // namespace sketcher
} // namespace schrodinger
