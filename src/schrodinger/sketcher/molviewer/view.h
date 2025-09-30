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
class SketcherModel;

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
    void setSketcherModel(SketcherModel* sketcher_model);

    /**
     * @return true if the user is in the middle of a pinch gesture rotation or
     * translation.  False otherwise.
     */
    bool isDuringPinchGesture();

  public slots:
    /**
     * Scale and center the view matrix so that objects in the scene are
     * visible, but avoid zooming in too much for small molecules
     * @param selection_only If true, only consider selected objects
     */
    void fitToScreen(bool selection_only = false);

    void fitAllToScreen()
    {
        fitToScreen(false);
    }

    /**
     * If needed scale the view matrix to include objects that are off-screen.
     */
    void zoomOutToIncludeAll();

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

    /**
     * Emitted when a pinch gesture has completed
     */
    void pinchGestureFinished();

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
     * Override QGraphicsView::showEvent so that we can call fitToScreen when
     * the view is first shown if necessary
     */
    void showEvent(QShowEvent* event) override;

    /**
     * Make sure that the scene rectangle contains everything that is currently
     * displayed in the viewport (including empty space). This prevents visual
     * artefacts when the selection rectangle or lasso is drawn outside of the
     * scene rect.
     */
    void enlargeSceneIfNeeded();

    /**
     * Center the view on the given point
     */
    void centerViewportOn(QPointF point);

    /**
     * Scale the matrix, multiply the current zoom level by scale_factor, but
     * cap the  result  at 1, so that the scene is never zoomed in more than
     * it's starting factor. All view zooming should be done using this function
     * to avoid zooming in too close.
     */
    void scaleSafely(qreal scale_factor);

    /**
     * Utility function to fit given rectangle to the screen.
     */
    void fitRecToScreen(const QRectF& rec);

    /**
     * Update the cursor when the color scheme is changed, since the cursor and
     * the background should always be contrasting colors
     */
    void onNewCursorColorRequested();

    /**
     * Update the cursor using the specified cursor hint
     */
    void updateCursor(const QPixmap& cursor_hint);

    MolModel* m_mol_model = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
    QPixmap m_cursor_hint;
    bool m_currently_pinching_trackpad = false;
    bool m_delayed_fit_to_screen = false;
    bool m_initial_geometry_set = false;
};

} // namespace sketcher
} // namespace schrodinger
