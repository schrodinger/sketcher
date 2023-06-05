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
    AbstractSceneTool(scene, mol_model)
{
}

void RotateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    // calculate angle increment of the event,
    // considering the centroid of the molecule as the origin
    auto centroid = to_scene_xy(m_mol_model->findCentroid());
    QLineF current_position_line(centroid, event->scenePos());
    QLineF last_position_line(centroid, event->lastScenePos());
    auto angle = last_position_line.angleTo(current_position_line);
    m_mol_model->rotateByAngle(angle);
}

} // namespace sketcher
} // namespace schrodinger