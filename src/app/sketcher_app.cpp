#include <QApplication>

#include "schrodinger/sketcher/sketcher_widget.h"

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
    schrodinger::sketcher::SketcherWidget sketcher_widget;
    sketcher_widget.show();
    return application.exec();
}
