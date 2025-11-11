#include "schrodinger/sketcher/molviewer/view.h"

#include <algorithm>

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

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

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

    /**
     * Set the viewport (since we don't have scrollbars, this is the same as
     * the scene rect for now)
     */
    setSceneRect(mapToScene(rect()).boundingRect());

    // We don't need to call enlargeSceneIfNeeded here since the View doesn't
    // have a size yet, so we're guaranteed to get a resizeEvent call before
    // View is painted.

    setBackgroundBrush(QBrush(LIGHT_BACKGROUND_COLOR));
}

View::View(QWidget* parent) : View(nullptr, parent)
{
}

bool View::isDuringPinchGesture()
{
    return m_currently_pinching_trackpad;
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
    enlargeSceneIfNeeded();
}

void View::pinchTriggered(QPinchGesture* gesture)
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene || !m_mol_model) {
        return;
    }

    auto state = gesture->state();
    if (state == Qt::GestureStarted) {
        m_currently_pinching_trackpad = true;
    } else {
        // zoom
        auto scale_factor = gesture->scaleFactor();
        scaleSafely(scale_factor);

        auto angle = gesture->rotationAngle() - gesture->lastRotationAngle();

        auto center_of_rotation = find_centroid(
            *(m_mol_model->getMol()), m_mol_model->getNonMolecularObjects());
        m_mol_model->rotateByAngle(-angle, center_of_rotation);
        if (state == Qt::GestureFinished) {
            m_currently_pinching_trackpad = false;
            emit pinchGestureFinished();
        }
    }
}

void View::enlargeSceneIfNeeded()
{
    QGraphicsScene* cur_scene = scene();
    if (!cur_scene) {
        return;
    }
    QRectF scene_rect = cur_scene->sceneRect();
    QRectF view_rect = mapToScene(rect()).boundingRect();
    cur_scene->setSceneRect(scene_rect.united(view_rect));
}

void View::centerViewportOn(QPointF point)
{
    auto scene_rect = sceneRect();
    scene_rect.moveCenter(point);
    setSceneRect(scene_rect);
}

void View::wheelEvent(QWheelEvent* event)
{
    // convert wheel movement to a scaling factor, using a function of type y =
    // K * 2^x, where K is an empirical constant to make the zooming feel
    // natural
    qreal scal = pow(2.0, event->angleDelta().y() / 2400.0);
    scaleSafely(scal);
}

void View::zoomOutToIncludeAll()
{
    Scene* cur_scene = dynamic_cast<Scene*>(scene());
    if (!cur_scene) {
        return;
    }
    auto visible_rect = mapToScene(rect()).boundingRect();
    QRectF rec = cur_scene->getSceneItemsBoundingRect();
    rec = rec.united(visible_rect);
    if (rec == visible_rect) {
        // If the bounding rect is the same as the visible rect, then we
        // don't need to zoom out, so just return
        return;
    }
    fitRecToScreen(rec);
}

void View::fitToScreen(bool selection_only)
{
    if (!m_initial_geometry_set) {
        // If the view hasn't been shown yet then it doesn't know how large it
        // will be, so this method would fit the molecule into the initial
        // default widget size, which is tiny. The zoom doesn't get updated when
        // the window is shown, so we'd wind up with a normal size window and a
        // very tiny molecule. To avoid this, we don't do anything now and
        // instead re-call this method as soon as the view is shown.
        m_delayed_fit_to_screen = true;
        return;
    }
    Scene* cur_scene = dynamic_cast<Scene*>(scene());
    if (!cur_scene) {
        return;
    }
    QRectF rec = cur_scene->getInteractiveItemsBoundingRect(
        InteractiveItemFlag::ALL, selection_only);
    // SKETCH-1703 SKETCH-2534 if the sketcher is editable, make the bounding
    // rect a bit bigger to avoid having the molecule too close to the border
    if (m_sketcher_model && !m_sketcher_model->isSelectOnlyModeActive()) {
        rec.adjust(-rec.width() * FIT_TO_SCREEN_MARGIN_FACTOR,
                   -rec.height() * FIT_TO_SCREEN_MARGIN_FACTOR,
                   rec.width() * FIT_TO_SCREEN_MARGIN_FACTOR,
                   rec.height() * FIT_TO_SCREEN_MARGIN_FACTOR);
    }
    fitRecToScreen(rec);
    // After fitting to screen, we want to ensure that the scene is not too
    // zoomed in
    scaleSafely(1.0);
}

void View::fitRecToScreen(const QRectF& rec)
{

    if (!rec.isValid()) {
        return;
    }
    fitInView(rec, Qt::KeepAspectRatio);
    centerViewportOn(rec.center());
    enlargeSceneIfNeeded();
}

void View::setMolModel(MolModel* mol_model)
{
    m_mol_model = mol_model;
    connect(mol_model, &MolModel::newMoleculeAdded, this,
            &View::fitAllToScreen);
    connect(mol_model, &MolModel::modelChanged, this,
            &View::zoomOutToIncludeAll);
}

void View::setSketcherModel(SketcherModel* sketcher_model)
{
    m_sketcher_model = sketcher_model;
    connect(sketcher_model, &SketcherModel::displaySettingsChanged, this,
            &View::onNewCursorColorRequested);
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
    enlargeSceneIfNeeded();
}

void View::keyPressEvent(QKeyEvent* event)
{
    auto distance = VIEW_SCALE * KEY_SCROLL_BOND_LENGTH_RATIO;
    switch (Qt::Key(event->key())) {
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
        default:
            break;
    }
    event->accept();
    QGraphicsView::keyPressEvent(event);
}

void View::leaveEvent(QEvent* event)
{
    QGraphicsView::leaveEvent(event);
    QGraphicsScene* cur_scene = scene();
    if (auto scene = dynamic_cast<Scene*>(cur_scene)) {
        scene->onMouseLeave();
    }
}

void View::showEvent(QShowEvent* event)
{
    QGraphicsView::showEvent(event);
    m_initial_geometry_set = true;
    if (m_delayed_fit_to_screen) {
        // This view got a fitToScreen call before it was shown (which means
        // that fitToScreen wouldn't have known how large the view was going to
        // be), so run fitToScreen now that geometry has been set
        fitToScreen();
        m_delayed_fit_to_screen = false;
    }
}

void View::onNewCursorColorRequested()
{
    updateCursor(m_cursor_hint);
}

void View::onNewCursorHintRequested(const QPixmap& cursor_hint)
{
    m_cursor_hint = cursor_hint;
    updateCursor(cursor_hint);
}

void View::updateCursor(const QPixmap& cursor_hint)
{
    // get the arrow image, colored for either light or dark mode
    QColor arrow_color, outline_color;
    if (m_sketcher_model->hasDarkColorScheme()) {
        arrow_color = DARK_CURSOR_COLOR;
        outline_color = DARK_BACKGROUND_COLOR;
    } else {
        arrow_color = LIGHT_CURSOR_COLOR;
        outline_color = LIGHT_BACKGROUND_COLOR;
    }
    auto arrow = get_arrow_cursor_pixmap(arrow_color, outline_color);

    // Figure out the size required to fit both the arrow and the cursor hint.
    // Note that, without the +1, the right-most column and bottom-most row of
    // inset's pixels will be cut off.
    int combined_width =
        std::max(arrow.width(), CURSOR_HINT_X + cursor_hint.width() + 1);
    int combined_height =
        std::max(arrow.height(), CURSOR_HINT_Y + cursor_hint.height() + 1);
    QPixmap combined = QPixmap(combined_width, combined_height);

    // paint both the arrow and the cursor hint into a single pixmap
    combined.fill(Qt::transparent);
    {
        QPainter painter(&combined);
        painter.drawPixmap(0, 0, arrow);
        painter.drawPixmap(CURSOR_HINT_X, CURSOR_HINT_Y, cursor_hint);
    } // destroy the painter to finish painting

    QCursor cursor(combined, CURSOR_HOTSPOT_X, CURSOR_HOTSPOT_Y);
    setCursor(cursor);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/view.moc"
