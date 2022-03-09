#include <QString>

namespace schrodinger
{
namespace sketcher
{

// TODO: These are LiveDesign colors; add Maestro colors in SKETCH-1342
const QString TOOL_BUTTON_STYLE{
    "QToolButton:checked { background-color: #b8d7ac; }\n"
    "QToolButton:hover:!checked { background-color: #daebd5 ; }\n"};

} // sketcher
} // schrodinger
