#include "schrodinger/sketcher/tool/arrow_plus_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger
{
namespace sketcher
{

ArrowPlusSceneTool::ArrowPlusSceneTool(const NonMolecularType& type,
                                       Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_type(type)
{
}

void ArrowPlusSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    auto coords = to_mol_xy(event->scenePos());
    m_mol_model->addNonMolecularObject(m_type, coords);
}

QPixmap ArrowPlusSceneTool::createDefaultCursorPixmap() const
{
    QString path = m_type == NonMolecularType::RXN_ARROW
                       ? ":/icons/reaction_arrow.svg"
                       : ":/icons/reaction_plus.svg";
    return cursor_hint_from_svg(path);
}

} // namespace sketcher
} // namespace schrodinger
