#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

AbstractSceneTool::AbstractSceneTool(Scene* scene, MolModel* mol_model) :
    m_scene(scene),
    m_mol_model(mol_model)
{
}

void AbstractSceneTool::onMousePress(QGraphicsSceneMouseEvent* const event)
{
    m_mouse_pressed = true;
    m_mouse_press_scene_pos = event->scenePos();
}

void AbstractSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMouseRelease(QGraphicsSceneMouseEvent* const event)
{
    m_mouse_pressed = false;
    m_mouse_press_scene_pos = QPointF();
}

void AbstractSceneTool::onMouseLeave()
{
}

void AbstractSceneTool::onDragStart(QGraphicsSceneMouseEvent* const event)
{
    m_drag_started = true;
}

void AbstractSceneTool::onDragMove(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onDragRelease(QGraphicsSceneMouseEvent* const event)
{
    m_drag_started = false;
}

void AbstractSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
}

void AbstractSceneTool::onMouseDoubleClick(
    QGraphicsSceneMouseEvent* const event)
{
}

std::vector<QGraphicsItem*> AbstractSceneTool::getGraphicsItems()
{
    return {};
}

} // namespace sketcher
} // namespace schrodinger
