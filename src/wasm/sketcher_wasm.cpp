#include <emscripten.h>
#include <emscripten/bind.h>

#include <QApplication>
#include <QFile>

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

void apply_stylesheet(QApplication& app)
{
    QFile styleFile(":resources/wasm_stylesheet.qss");
    styleFile.open(QFile::ReadOnly);
    QString style(styleFile.readAll());
    app.setStyleSheet(style);
}

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
    apply_stylesheet(application);

    auto& sk = get_sketcher_instance();
    sk.show();
    return application.exec();
}
