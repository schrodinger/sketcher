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

void TranslateSceneTool::onDragStart(QGraphicsSceneMouseEvent* event)
{
    // if the translate drag starts inside a selection, move only the selction
    if (m_scene->getSelectionRect().contains(m_mouse_press_scene_pos)) {
        setCurrentSelectionAsObjectsToMove();
    }
    MoveRotateSceneTool::onDragStart(event);
}

void TranslateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    translate(event, m_atoms_to_move, m_non_mol_objs_to_move);
}

void TranslateSceneTool::onMouseClick(QGraphicsSceneMouseEvent* event)
{
    /*
     * clicking with the right mouse button shows the context menu
     */
    m_scene->showContextMenu(event);
}

} // namespace sketcher
} // namespace schrodinger