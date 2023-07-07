#include "schrodinger/sketcher/tool/arrow_plus_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

ArrowPlusSceneTool::ArrowPlusSceneTool(const NonMolecularType& type,
                                       Scene* scene, MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model),
    m_type(type)
{
}

void ArrowPlusSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    auto coords = to_mol_xy(event->scenePos());
    m_mol_model->addNonMolecularObject(m_type, coords);
}

} // namespace sketcher
} // namespace schrodinger
