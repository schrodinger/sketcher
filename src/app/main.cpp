// @copyright Schrodinger, LLC - All Rights Reserved

#include <QApplication>
#include <QtPlugin>

#include "schrodinger/sketcher/example_widget.h"

// Q_IMPORT_PLUGIN(QXcbIntegrationPlugin)

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
    schrodinger::sketcher::ExampleWidget widget;
    widget.show();
    return application.exec();
}
