#pragma once

#include <vector>

#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPointF>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;

/**
 * A base class for scene tools, which encapsulate mouse behavior in the Scene
 * for a specific Sketcher tool (e.g. the draw nitrogen tool, or the erase
 * tool).  When the user changes their tool in Sketcher, the Scene will
 * instantiate a new scene tool, as well as add all of the graphics items
 * returned from getGraphicsItems().  These graphics items will be removed from
 * the Scene when the user changes to a new tool, so destroying the graphics
 * items is the responsibility of the scene tool.
 *
 * When the user presses, moves, or releases the mouse, the corresponding method
 * will be called with the associated event.  Same with mouse drags and clicks.
 * A mouse drag starts when the user has moved the mouse
 * QApplication::startDragDistance() pixels.  A mouse click is any press and
 * release in which no drag was started.
 *
 * Subclasses should reimplement some or all of the virtual methods.  Overridden
 * methods should always start by calling the super-class method to ensure that
 * bookkeeping variables (m_mouse_pressed, m_drag_started, and
 * m_mouse_press_scene_pos) are correctly updated.
 */
class SKETCHER_API AbstractSceneTool
{
  public:
    AbstractSceneTool(Scene* scene, MolModel* mol_model);
    virtual ~AbstractSceneTool() = default;

    /**
     * Called whenever the user presses the left mouse button
     */
    virtual void onMousePress(QGraphicsSceneMouseEvent* const event);

    /**
     * Called whenever the user moves the mouse
     */
    virtual void onMouseMove(QGraphicsSceneMouseEvent* const event);

    /**
     * Called whenever the user releases the left mouse button
     */
    virtual void onMouseRelease(QGraphicsSceneMouseEvent* const event);

    /**
     * Called whenever the user starts a drag by moving the mouse
     * QApplication::startDragDistance() pixels with the left mouse button held
     * down.  Note that onDragMove() will also be called with the same mouse
     * event.
     */
    virtual void onDragStart(QGraphicsSceneMouseEvent* const event);

    /**
     * Called when the mouse is moved during a click-and-drag
     */
    virtual void onDragMove(QGraphicsSceneMouseEvent* const event);

    /**
     * Called when the user releases the mouse button, ending a click-and-drag
     */
    virtual void onDragRelease(QGraphicsSceneMouseEvent* const event);

    /**
     * Called when the user clicks the left mouse button.  A mouse click is any
     * press and release in which no drag was started.
     */
    virtual void onMouseClick(QGraphicsSceneMouseEvent* const event);

    /**
     * @return all graphics items that should be added to the scene while this
     * tool is active (e.g. predictive highlighting, marquee selection outline).
     * These graphics items will be removed from the Scene when the user changes
     * to a new tool, so destroying the graphics items is the responsibility of
     * this class.
     */
    virtual std::vector<QGraphicsItem*> getGraphicsItems();

    // TODO
    // QCursor getCursor();

  protected:
    Scene* m_scene;
    MolModel* m_mol_model;
    bool m_mouse_pressed = false;
    bool m_drag_started = false;
    QPointF m_mouse_press_scene_pos = QPointF();
};

} // namespace sketcher
} // namespace schrodinger
