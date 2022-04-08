#include <QString>

namespace schrodinger
{
namespace sketcher
{

const QString DISABLED_COLOR{"#b6b6b6"};

const QString TOOL_BUTTON_WITH_POPUP_BASE_STYLE{"QToolButton {"
                                                " border : none;"
                                                " padding: 0px;"
                                                " margin: 0px;"
                                                " spacing: 0px;"
                                                "}"};

const QString TOOL_BUTTON_WITH_POPUP_AND_ARROW_STYLE =
    TOOL_BUTTON_WITH_POPUP_BASE_STYLE +
    QString{"QToolButton {"
            " background-image:url(\":icons/corner_arrow.png\"); "
            " background-position: bottom right; "
            " background-repeat: no-repeat; "
            "}"
            "QToolButton:disabled {"
            " background-image:url(\":icons/disabled_corner_arrow.png\"); "
            "}"};

// TODO: These are LiveDesign colors; add Maestro colors in SKETCH-1342
const QString TOOL_BUTTON_STYLE{
    "QToolButton:checked { background-color: #b8d7ac; }\n"
    "QToolButton:hover:!checked { background-color: #daebd5 ; }\n"};

// TODO: Extract the general Sketcher background color into a constant.
const QString WINDOW_STYLE{"QMainWindow { background-color: white; }\n"};
} // namespace sketcher
} // namespace schrodinger
