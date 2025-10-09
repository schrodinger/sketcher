#pragma once

#include <vector>

#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QObject>
#include <QPointF>

#include "schrodinger/sketcher/definitions.h"

class QPixmap;

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;

/**
 * A base class for scene tools, which encapsulate mouse behavior in the Scene
 * for a specific Sketcher tool (e.g. the draw nitrogen tool, or the erase
 * tool). When the user changes their tool in Sketcher, the Scene will
 * instantiate a new scene tool, as well as add all of the graphics items
 * returned from getGraphicsItems(). These graphics items will be removed from
 * the Scene when the user changes to a new tool, so destroying the graphics
 * items is the responsibility of the scene tool.
 *
 * When the user presses, moves, or releases the mouse, the corresponding method
 * will be called with the associated event. Same with mouse drags and clicks.
 * A mouse drag starts when the user has moved the mouse
 * QApplication::startDragDistance() pixels. A mouse click is any press and
 * release in which no drag was started.
 *
 * Subclasses should reimplement some or all of the virtual methods. Overridden
 * methods should always start by calling the super-class method to ensure that
 * bookkeeping variables (m_mouse_pressed, m_drag_started, and
 * m_mouse_press_scene_pos) are correctly updated.
 */
class SKETCHER_API AbstractSceneTool : public QObject
{
    Q_OBJECT
  public:
    AbstractSceneTool(Scene* scene, MolModel* mol_model);

    /**
     * Event handlers for mouse events. Note that these event handlers are
     * called by the scene, and are *not* called directly by Qt. (This class
     * isn't a QWidget and doesn't have a parent, so Qt will never send it any
     * events directly.) To change mouse behavior, subclasses should override
     * the onWhateverButtonSomethinged methods below, not these methods.
     */
    void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event);

    /**
     * @return all graphics items that should be added to the scene while this
     * tool is active (e.g. predictive highlighting, marquee selection outline).
     * These graphics items will be removed from the Scene when the user changes
     * to a new tool, so destroying the graphics items is the responsibility of
     * this class.
     */
    virtual std::vector<QGraphicsItem*> getGraphicsItems();

    /**
     * @return a cached pixmap of a hint image to attach to the mouse cursor.
     * Subclasses should override createDefaultCursorPixmap to provide a custom
     * hint image.
     */
    QPixmap getDefaultCursorPixmap() const;

    /**
     * Called when the structure of the molecule changes. It does nothing
     * by default, but can be reimplemented in subclasses to update the tool's
     * state (e.g. graphics items that needs to be updated)
     */
    virtual void onStructureUpdated();

    /**
     * Called when the background color of the scene changes. It does nothing
     * by default, but can be reimplemented in subclasses to update the tool's
     * state (e.g. graphics items that are on different colors on dark
     * background)
     */
    virtual void updateColorsAfterBackgroundColorChange(bool is_dark_mode);

    /**
     * Called when the selection of the scene changes. It does nothing
     * by default, but can be reimplemented in subclasses to update the tool's
     * state (e.g. graphics items that needs to be updated)
     */
    virtual void onSelectionChanged();

    /**
     * Called whenever the mouse leaves the scene.
     */
    virtual void onMouseLeave();

  signals:
    /**
     * A signal emitted when the context menu should be shown.
     * @param event the mouse event that triggered the context menu request
     */
    void contextMenuRequested(QGraphicsSceneMouseEvent* event);

    /**
     * A signal emitted when the cursor hint should be changed
     * @param pixmap the new cursor hint pixmap
     */
    void newCursorHintRequested(const QPixmap& pixmap);

    /**
     * A signal emitted when a rotation or translation drag is started.
     */
    void atomDragStarted();

    /**
     * A signal emitted when a rotation or translation drag is finished.
     * @param were_atoms_merged whether or not the drag resulted in atoms being
     * merged
     */
    void atomDragFinished(const bool were_atoms_merged);

  protected:
    /**
     * @return a pixmap of a hint image to attach to the mouse cursor.
     * Subclasses should override this method to provide a custom hint image.
     * The default implementation returns a null pixmap, which clears the cursor
     * hint.
     */
    virtual QPixmap createDefaultCursorPixmap() const;

    /**
     * Called whenever the user moves the mouse
     */
    virtual void onMouseMove(QGraphicsSceneMouseEvent* const event);

    /**
     * Called whenever the user presses or releases a mouse button. For behavior
     * that triggers on mouse clicks, override the onWhateverButtonClick methods
     * below instead of these methods.
     */
    virtual void onLeftButtonPress(QGraphicsSceneMouseEvent* const event);
    virtual void onLeftButtonRelease(QGraphicsSceneMouseEvent* const event);

    virtual void onMiddleButtonPress(QGraphicsSceneMouseEvent* const event);
    virtual void onMiddleButtonRelease(QGraphicsSceneMouseEvent* const event);

    virtual void onRightButtonPress(QGraphicsSceneMouseEvent* const event);
    virtual void onRightButtonRelease(QGraphicsSceneMouseEvent* const event);

    /**
     * Called whenever the user drags the mouse, which we define as moving the
     * mouse QApplication::startDragDistance() pixels with the mouse button held
     * down. Note that, for the mouse event that starts the drag, *both* the
     * onWhateverButtonDragStart() and onWhateverButtonDragMove() methods will
     * be called.
     */
    virtual void onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event);
    virtual void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event);
    virtual void onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event);

    virtual void onMiddleButtonDragStart(QGraphicsSceneMouseEvent* const event);
    virtual void onMiddleButtonDragMove(QGraphicsSceneMouseEvent* const event);
    virtual void
    onMiddleButtonDragRelease(QGraphicsSceneMouseEvent* const event);

    virtual void onRightButtonDragStart(QGraphicsSceneMouseEvent* const event);
    virtual void onRightButtonDragMove(QGraphicsSceneMouseEvent* const event);
    virtual void
    onRightButtonDragRelease(QGraphicsSceneMouseEvent* const event);

    /**
     * Called when the user clicks a mouse button. A mouse click is any press
     * and release in which no drag was started.
     */
    virtual void onLeftButtonClick(QGraphicsSceneMouseEvent* const event);
    virtual void onMiddleButtonClick(QGraphicsSceneMouseEvent* const event);
    virtual void onRightButtonClick(QGraphicsSceneMouseEvent* const event);

    /**
     * Called when the user double clicks the left mouse button.
     */
    virtual void onLeftButtonDoubleClick(QGraphicsSceneMouseEvent* const event);

    Scene* m_scene;
    MolModel* m_mol_model;
    // Note that Qt::NoButton is falsey and LeftButton, MiddleButton, and
    // RightButton are all truthy. We use the precise button identity during
    // mouse drags in case the user presses multiple mouse buttons during the
    // drag.
    Qt::MouseButton m_mouse_pressed = Qt::NoButton;
    bool m_drag_started = false;
    bool m_is_during_double_click = false;
    QPointF m_mouse_press_scene_pos = QPointF();
    QPointF m_drag_start_scene_pos = QPointF();
    // the cursor hint is created during the first getDefaultCursorPixmap call,
    // so we mark these values as mutable since that method is const
    mutable QPixmap m_default_cursor_hint;
    mutable bool m_created_default_cursor_hint = false;

  private:
    // the screen pos (as opposed to the *scene* pos) is private since
    // subclasses are far more likely to access it accidentally as a typo of
    // m_mouse_press_scene_pos than they are to have a legitimate use for the
    // screen pos. This base class uses the screen pos (i.e. this value) to
    // decide when a drag has started.
    QPointF m_mouse_press_screen_pos = QPointF();
};

} // namespace sketcher
} // namespace schrodinger
