#include <QApplication>

#include "schrodinger/sketcher/example_widget.h"

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
    schrodinger::sketcher::ExampleWidget widget;
    widget.show();
    return application.exec();
}
