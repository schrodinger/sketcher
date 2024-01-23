#include "schrodinger/sketcher/tool/rotate_scene_tool.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{
RotateSceneTool::RotateSceneTool(Scene* scene, MolModel* mol_model) :
    MoveRotateSceneTool(scene, mol_model)
{
}

void RotateSceneTool::onDragStart(QGraphicsSceneMouseEvent* event)
{
    // if a selection is present, rotate only the selection
    setCurrentSelectionAsObjectsToMove();
    MoveRotateSceneTool::onDragStart(event);
}

void RotateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    auto center_of_rotation = find_centroid(
        *(m_mol_model->getMol()), m_atoms_to_move, m_non_mol_objs_to_move);
    rotate(event, to_scene_xy(center_of_rotation), m_atoms_to_move,
           m_non_mol_objs_to_move);
}

} // namespace sketcher
} // namespace schrodinger