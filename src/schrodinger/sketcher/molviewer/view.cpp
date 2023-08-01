#include "schrodinger/sketcher/molviewer/view.h"

#include <QGraphicsScene>
#include <QPainter>
#include <QResizeEvent>
#include <QPinchGesture>
#include <QGestureEvent>
#include <QPinchGesture>
#include <QGestureEvent>
#include <Qt>
#include <QtGlobal>
#include <QWheelEvent>
#include <QWheelEvent>

#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

View::View(QGraphicsScene* scene, QWidget* parent) :
    QGraphicsView(scene, parent)
{
    setMouseTracking(true);

    grabGesture(Qt::PinchGesture);
    setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);

    // disable scrollbars
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

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

bool View::event(QEvent* event)
{
    if (event->type() == QEvent::Gesture) {
        return gestureEvent(static_cast<QGestureEvent*>(event));
    }
    return QGraphicsView::event(event);
}

bool View::gestureEvent(QGestureEvent* event)
{
    if (QGesture* pinch = event->gesture(Qt::PinchGesture)) {
        pinchTriggered(static_cast<QPinchGesture*>(pinch));
        event->accept();
    }
    return true;
}

void View::scaleSafely(qreal scale_factor)
{
    scale(scale_factor, scale_factor);
    float zoom_threshold = 1.0;
    auto matrix = transform();
    float m11 = matrix.m11();
    float m22 = matrix.m22();
    if (m11 > zoom_threshold || m22 > zoom_threshold) {
        matrix.setMatrix(zoom_threshold, matrix.m12(), matrix.m13(),
                         matrix.m21(), zoom_threshold, matrix.m23(),
                         matrix.m31(), matrix.m32(), matrix.m33());
        setTransform(matrix);
    }
}

void View::pinchTriggered(QPinchGesture* gesture)
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene || !m_mol_model) {
        return;
    }

    auto state = gesture->state();
    if (state != Qt::GestureStarted) {
        // zoom
        auto scale_factor = gesture->scaleFactor();
        scaleSafely(scale_factor);

        auto angle = gesture->rotationAngle() - gesture->lastRotationAngle();

        auto center_of_rotation = find_centroid(
            *(m_mol_model->getMol()), m_mol_model->getNonMolecularObjects());
        m_mol_model->rotateByAngle(-angle, center_of_rotation);
    }
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

void View::adjustSceneAroundItems()
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene) {
        return;
    }
    // get the bounding rectangle of all items in the scene
    QRectF items_bounding_rect = cur_scene->itemsBoundingRect();
    if (!items_bounding_rect.isValid()) {
        return;
    }
    cur_scene->setSceneRect(items_bounding_rect);
}

void View::wheelEvent(QWheelEvent* event)
{
    // convert wheel movement to a scaling factor, using a function of type y =
    // K * 2^x, where K is an empirical constant to make the zooming feel
    // natural
    qreal scal = pow(2.0, event->angleDelta().y() / 2400.0);
    scaleSafely(scal);
}

void View::fitToScreen()
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene) {
        return;
    }

    QRectF rec = cur_scene->itemsBoundingRect();
    // SKETCH-1703 make the bounding rect a bit bigger to avoid having the
    // molecule too close to the border
    rec.adjust(-rec.width() * FIT_TO_SCREEN_MARGIN_FACTOR,
               -rec.height() * FIT_TO_SCREEN_MARGIN_FACTOR,
               rec.width() * FIT_TO_SCREEN_MARGIN_FACTOR,
               rec.height() * FIT_TO_SCREEN_MARGIN_FACTOR);

    if (rec.isValid()) {
        fitInView(rec, Qt::KeepAspectRatio);
        // avoid zooming in too much
        scaleSafely(1.0);
    }
    adjustSceneAroundItems();
    enlargeSceneIfNeeded();
}

void View::setMolModel(schrodinger::sketcher::MolModel* mol_model)
{
    m_mol_model = mol_model;
}

void View::translateViewportFromScreenCoords(
    const QPointF& start_screen_position, const QPointF& end_screen_position)
{
    QPointF start_position = mapToScene(start_screen_position.toPoint());
    QPointF end_position = mapToScene(end_screen_position.toPoint());
    translateViewport(start_position - end_position);
}

void View::translateViewport(const QPointF& delta)
{
    setSceneRect(sceneRect().translated(delta));
}

void View::keyPressEvent(QKeyEvent* event)
{
    auto distance = VIEW_SCALE * KEY_SCROLL_BOND_LENGTH_RATIO;
    switch (event->key()) {
        case Qt::Key_Up:
            translateViewport(QPointF(0, distance));
            break;
        case Qt::Key_Down:
            translateViewport(QPointF(0, -distance));
            break;
        case Qt::Key_Right:
            translateViewport(QPointF(-distance, 0));
            break;
        case Qt::Key_Left:
            translateViewport(QPointF(distance, 0));
            break;
    }
    event->accept();
}

void View::leaveEvent(QEvent* event)
{
    QGraphicsView::leaveEvent(event);
    QGraphicsScene* cur_scene = scene();
    if (auto scene = dynamic_cast<Scene*>(cur_scene)) {
        scene->onMouseLeave();
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/view.moc"
