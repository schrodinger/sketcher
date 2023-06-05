#include "schrodinger/sketcher/tool/translate_scene_tool.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{
TranslateSceneTool::TranslateSceneTool(Scene* scene, MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model)
{
}

void TranslateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    auto distance = event->scenePos() - event->lastScenePos();
    m_mol_model->translateByVector(to_mol_xy(distance));
}

void TranslateSceneTool::onMouseClick(QGraphicsSceneMouseEvent* event)
{
    m_scene->showContextMenu(event);
}

} // namespace sketcher
} // namespace schrodinger