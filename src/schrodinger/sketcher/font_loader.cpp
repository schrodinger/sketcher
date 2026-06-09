#include "schrodinger/sketcher/font_loader.h"

#include <fmt/format.h>

#include <QDebug>
#include <QFile>
#include <QFontDatabase>

namespace schrodinger
{
namespace sketcher
{

void load_font_resources()
{
    static bool already_loaded = false;
    if (already_loaded) {
        // this function only needs to run once per process, so we can return
        // early if this isn't the first time it's been called
        return;
    }
    for (const auto& font_path : {":/resources/fonts/Arimo-Regular.ttf",
                                  ":/resources/fonts/Arimo-Bold.ttf",
                                  ":/resources/fonts/Arimo-Italic.ttf",
                                  ":/resources/fonts/Arimo-BoldItalic.ttf"}) {
        if (!QFile::exists(font_path)) {
            // The font resource isn't embedded in the binary, which almost
            // always means a typo in the path above or a problem compiling the
            // .qrc file into the application. This is a programmer error, so
            // fail loudly.
            throw std::runtime_error(
                fmt::format("Font resource not found: {}", font_path));
        }
        if (QFontDatabase::addApplicationFont(font_path) == -1) {
            // The font data is present but the OS refused to register it. This
            // happens on machines (commonly Windows enterprise installations)
            // where a security policy blocks loading fonts that aren't
            // installed in the protected system fonts directory. In this
            // scenario, warn and fall back to a system font.
            qWarning() << "Failed to register font (resource data present):"
                       << font_path;
        }
    }
    already_loaded = true;
}

} // namespace sketcher
} // namespace schrodinger
