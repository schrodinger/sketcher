#include "schrodinger/sketcher/font_loader.h"

#include <fmt/format.h>

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
        if (QFontDatabase::addApplicationFont(font_path) == -1) {
            throw std::runtime_error(
                fmt::format("Failed to load font: {}", font_path));
        }
    }
    already_loaded = true;
}

} // namespace sketcher
} // namespace schrodinger
