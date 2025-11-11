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

/**
 * Whether the Sketcher can be used to draw atomistic models (i.e. small
 * molecules), monomeric models (i.e. HELM models), or both
 */
typedef int InterfaceTypeType;
namespace InterfaceType
{
enum InterfaceTypeType {
    ATOMISTIC = 1 << 0,
    MONOMERIC = 1 << 1,
    ATOMISTIC_OR_MONOMERIC = ATOMISTIC | MONOMERIC,
};
}

} // namespace sketcher
} // namespace schrodinger
