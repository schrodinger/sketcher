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
    "QToolButton:checked { background-color: rgba(184, 215, 172, .65); }\n"
    "QToolButton:hover:!checked { background-color: rgba(218, 235, 213, .5) ; "
    "}\n"};

// TODO: Extract the general Sketcher background color into a constant.
const QString GENERAL_STYLE{
    "QMainWindow { background-color: white; }\n"
    "schrodinger--sketcher--ModularPopup { background-color: white }"};
} // namespace sketcher
} // namespace schrodinger
