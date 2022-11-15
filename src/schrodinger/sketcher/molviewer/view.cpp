#include "schrodinger/sketcher/molviewer/view.h"

#include <Qt>
#include <QGraphicsScene>
#include <QPainter>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

View::View(QGraphicsScene* scene, QWidget* parent) :
    QGraphicsView(scene, parent)
{
    setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
}

View::View(QWidget* parent) : QGraphicsView(parent)
{
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/view.moc"
