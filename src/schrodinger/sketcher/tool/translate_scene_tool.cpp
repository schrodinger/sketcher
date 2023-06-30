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
    MoveRotateSceneTool(scene, mol_model)
{
}

void TranslateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    translate(event);
}

void TranslateSceneTool::onMouseClick(QGraphicsSceneMouseEvent* event)
{
    m_scene->showContextMenu(event);
}

} // namespace sketcher
} // namespace schrodinger