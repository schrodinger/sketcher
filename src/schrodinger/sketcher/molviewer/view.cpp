#include "schrodinger/sketcher/molviewer/view.h"

#include <Qt>
#include <QtGlobal>
#include <QGraphicsScene>
#include <QPainter>
#include <QResizeEvent>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

View::View(QGraphicsScene* scene, QWidget* parent) :
    QGraphicsView(scene, parent)
{
    setMouseTracking(true);
    setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);

    // We don't need to call enlargeSceneIfNeeded here since the View doesn't
    // have a size yet, so we're guaranteed to get a resizeEvent call before
    // View is painted.
}

View::View(QWidget* parent) : View(nullptr, parent)
{
}

void View::resizeEvent(QResizeEvent* event)
{
    QGraphicsView::resizeEvent(event);
    enlargeSceneIfNeeded();
    emit resized();
}

void View::enlargeSceneIfNeeded()
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene) {
        return;
    }
    QRectF scene_rect = cur_scene->sceneRect();
    QSizeF view_size = size();
    qreal extra_width = view_size.width() - scene_rect.width();
    extra_width = qMax(0.0, extra_width);
    qreal extra_height = view_size.height() - scene_rect.height();
    extra_height = qMax(0.0, extra_height);
    if (extra_width == 0.0 && extra_height == 0.0) {
        // The scene is already at least as large as the view
        return;
    }
    qreal half_width = extra_width / 2;
    qreal half_height = extra_height / 2;
    scene_rect.adjust(-half_width, -half_height, half_width, half_height);
    cur_scene->setSceneRect(scene_rect);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/view.moc"
