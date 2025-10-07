#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

#include <QApplication>
#include <QPixmap>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

AbstractSceneTool::AbstractSceneTool(Scene* scene, MolModel* mol_model) :
    QObject(),
    m_scene(scene),
    m_mol_model(mol_model)
{
}

/**
 * Call the appropriate scene tool member function for the given mouse button
 * @param scene_tool The scene tool instance to call the member function on
 * @param button The mouse button to use for deciding which method to call
 * @param event The Qt event that should be passed to the member function
 * @param left_button_method The member function to call if button is
 * Qt::LeftButton.
 * @param middle_button_method The member function to call if button is
 * Qt::MiddleButton.
 * @param right_button_method The member function to call if button is
 * Qt::RightButton.
 */
static void call_method_for_button(
    AbstractSceneTool* scene_tool, const Qt::MouseButton button,
    QGraphicsSceneMouseEvent* const event,
    void (AbstractSceneTool::*left_button_method)(QGraphicsSceneMouseEvent*),
    void (AbstractSceneTool::*middle_button_method)(QGraphicsSceneMouseEvent*),
    void (AbstractSceneTool::*right_button_method)(QGraphicsSceneMouseEvent*))
{
    if (button == Qt::LeftButton) {
        (scene_tool->*left_button_method)(event);
    } else if (button == Qt::MiddleButton) {
        (scene_tool->*middle_button_method)(event);
    } else if (button == Qt::RightButton) {
        (scene_tool->*right_button_method)(event);
    };
}

void AbstractSceneTool::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    if (m_mouse_pressed) {
        // there's already one mouse button pressed, so ignore any additional
        // buttons until the first one is released
        return;
    }
    m_mouse_press_scene_pos = event->scenePos();
    m_mouse_press_screen_pos = event->screenPos();
    auto btn = event->button();
    m_mouse_pressed = btn;
    call_method_for_button(this, btn, event,
                           &AbstractSceneTool::onLeftButtonPress,
                           &AbstractSceneTool::onMiddleButtonPress,
                           &AbstractSceneTool::onRightButtonPress);
}

void AbstractSceneTool::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // we used the m_mouse_pressed button here instead of checking
    // event->buttons() just in case the user decided to hold down two buttons
    // at once. This way, we'll always attribute the movement to the first
    // button that was pressed and ignore the second one.
    if (m_mouse_pressed != Qt::NoButton && !m_drag_started) {
        // check to see whether the mouse has moved far enough to start a
        // drag
        int drag_dist =
            (m_mouse_press_screen_pos - event->screenPos()).manhattanLength();
        if (drag_dist >= QApplication::startDragDistance()) {
            m_drag_started = true;
            m_drag_start_scene_pos = event->scenePos();
            call_method_for_button(this, m_mouse_pressed, event,
                                   &AbstractSceneTool::onLeftButtonDragStart,
                                   &AbstractSceneTool::onMiddleButtonDragStart,
                                   &AbstractSceneTool::onRightButtonDragStart);
        }
    }
    onMouseMove(event);
    if (m_drag_started) {
        call_method_for_button(this, m_mouse_pressed, event,
                               &AbstractSceneTool::onLeftButtonDragMove,
                               &AbstractSceneTool::onMiddleButtonDragMove,
                               &AbstractSceneTool::onRightButtonDragMove);
    }
}

void AbstractSceneTool::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    auto btn = event->button();
    if (m_is_during_double_click && btn == Qt::LeftButton) {
        // This event is the mouse release for the second click of a double
        // click, so we never got a corresponding mousePressEvent (since it
        // was a mouseDoubleClickEvent instead) and we've already handled
        // the double-click (in mouseDoubleClickEvent).
        m_is_during_double_click = false;
        return;
    }
    if (btn != m_mouse_pressed) {
        // if the user was holding down multiple buttons, we only care about the
        // first one that was pressed
        return;
    }
    if (m_drag_started) {
        call_method_for_button(this, btn, event,
                               &AbstractSceneTool::onLeftButtonDragRelease,
                               &AbstractSceneTool::onMiddleButtonDragRelease,
                               &AbstractSceneTool::onRightButtonDragRelease);
    } else {
        call_method_for_button(this, btn, event,
                               &AbstractSceneTool::onLeftButtonClick,
                               &AbstractSceneTool::onMiddleButtonClick,
                               &AbstractSceneTool::onRightButtonClick);
    }
    call_method_for_button(this, btn, event,
                           &AbstractSceneTool::onLeftButtonRelease,
                           &AbstractSceneTool::onMiddleButtonRelease,
                           &AbstractSceneTool::onRightButtonRelease);
    m_mouse_press_scene_pos = QPointF();
    m_mouse_press_screen_pos = QPointF();
    m_drag_started = false;
    m_mouse_pressed = Qt::NoButton;
}

void AbstractSceneTool::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event)
{
    // the !m_mouse_pressed check here ensures that we ignore double clicks that
    // happen while another mouse button is being held down
    if (event->button() == Qt::LeftButton && !m_mouse_pressed) {
        m_is_during_double_click = true;
        onLeftButtonDoubleClick(event);
    }
}

std::vector<QGraphicsItem*> AbstractSceneTool::getGraphicsItems()
{
    return {};
}

QPixmap AbstractSceneTool::createDefaultCursorPixmap() const
{
    return QPixmap();
}

QPixmap AbstractSceneTool::getDefaultCursorPixmap() const
{
    if (!m_created_default_cursor_hint) {
        m_default_cursor_hint = createDefaultCursorPixmap();
        m_created_default_cursor_hint = true;
    }
    return m_default_cursor_hint;
}

void AbstractSceneTool::onStructureUpdated()
{
}

void AbstractSceneTool::updateColorsAfterBackgroundColorChange(bool)
{
}

void AbstractSceneTool::onSelectionChanged()
{
}

void AbstractSceneTool::onMouseLeave()
{
}

void AbstractSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonPress(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonPress(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonPress(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonClick(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMiddleButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onRightButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onLeftButtonDoubleClick(
    QGraphicsSceneMouseEvent* const event)
{
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/tool/abstract_scene_tool.moc"
