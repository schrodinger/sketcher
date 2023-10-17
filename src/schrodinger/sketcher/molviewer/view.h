#pragma once

#include <QGraphicsView>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"

class QGraphicsScene;
class QResizeEvent;
class QGestureEvent;
class QPinchGesture;
class QWidget;

namespace schrodinger
{
namespace sketcher
{
class MolModel;

/**
 * A Qt graphics view for displaying molecules in a molviewer Scene.
 */
class SKETCHER_API View : public QGraphicsView
{
    Q_OBJECT
  public:
    View(QGraphicsScene* scene, QWidget* parent = nullptr);
    View(QWidget* parent = nullptr);
    void setMolModel(MolModel* mol_model);

  public slots:
    /**
     * Scale and center the view matrix so that all objects in the scene are
     * visible, but avoid zooming in too much for small molecules
     */
    void fitToScreen();

    /**
     * Translate the viewport by delta. This is used to pan the view around the
     * scene when the user hits the arrow keys.
     */
    void translateViewport(const QPointF& delta);

    /**
     * translate the viewport by the given amount in screen coordinates. This is
     * used to pan the view around the scene when the user uses the mouse tools
     */
    void translateViewportFromScreenCoords(const QPointF& start_screen_position,
                                           const QPointF& end_screen_position);

    /**
     * Update the cursor when the Scene requests a new cursor hint
     * @param cursor_hint A pixmap showing the desired cursor hint
     */
    void onNewCursorHintRequested(const QPixmap& cursor_hint);

  signals:
    void resized();

  protected:
    // Override the QGraphicsView method so we can call enlargeSceneIfNeeded
    void resizeEvent(QResizeEvent* event) override;

    // Override the QWidget method so we can detect trackpad gestures
    bool event(QEvent* event) override;

    // respond to a trackpad gesture, such as pinch to zoom and rotate
    bool gestureEvent(QGestureEvent* event);

    /** respond to a pinch gesture. This is called by gestureEvent and is
     * responsible for tracking the distance and angle of the user's fingers and
     * converting them into a zoom of the matrix and a rotation of the
     * molecule's coordinates respectively.
     */
    void pinchTriggered(QPinchGesture* gesture);

    // handle mouse wheel events, zooming in on wheel up and out on wheel down.
    void wheelEvent(QWheelEvent* event) override;

    // Override the QGraphicsView keyboard event method
    void keyPressEvent(QKeyEvent* event) override;

    // Override the QWidget method so that we can pass leave events to the Scene
    // (which passes them on to the scene tool)
    void leaveEvent(QEvent* event) override;

    /**
     * Make sure that the scene is at least as large as the view.  If the
     * scene is smaller than the view, then the view will center the scene.
     * Then, if a new item gets added that causes the scene to grow, the
     * view will re-center, which causes the molecule to jump around.
     */
    void enlargeSceneIfNeeded();

    /**
     * Make sure that the scene has enough space around the items so that it can
     * be centered.
     */
    void adjustSceneAroundItems();

    /**
     * Scale the matrix, multiply the current zoom level by scale_factor, but
     * cap the  result  at 1, so that the scene is never zoomed in more than
     * it's starting factor. All view zooming should be done using this function
     * to avoid zooming in too close.
     */
    void scaleSafely(qreal scale_factor);

    MolModel* m_mol_model = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
