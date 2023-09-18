#include <emscripten.h>
#include <emscripten/bind.h>

#include <QApplication>

#include "schrodinger/sketcher/sketcher_widget.h"

// For the WebAssembly build, we need to be able to get the sketcher
// instance we are running from a function/static method. We'll use a
// singleton for that.
WASMSketcher& get_sketcher_instance()
{
    static schrodinger::sketcher::SketcherWidget instance;
    return instance;
}

bool example_function()
{
    return true;
}

EMSCRIPTEN_BINDINGS(sketcher)
{
    emscripten::function("example_function", &example_function);
}

int main(int argc, char** argv)
{
    QApplication application(argc, argv);

    auto& sk = get_sketcher_instance();
    sk.show();
    return application.exec();
}
