#include "schrodinger/sketcher/tool/scene_tool_with_predictive_highlighting.h"

#include <QPointF>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

SceneToolWithPredictiveHighlighting::SceneToolWithPredictiveHighlighting(
    Scene* scene, MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model)
{
}

std::vector<QGraphicsItem*>
SceneToolWithPredictiveHighlighting::getGraphicsItems()
{
    return {&m_predictive_highlighting_item};
}

void SceneToolWithPredictiveHighlighting::onMouseMove(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractSceneTool::onMouseMove(event);
    if (!m_mouse_pressed) {
        QPointF scene_pos = event->scenePos();
        QGraphicsItem* item = m_scene->getTopInteractiveItemAt(scene_pos);
        m_predictive_highlighting_item.highlightItem(item);
    }
}

} // namespace sketcher
} // namespace schrodinger