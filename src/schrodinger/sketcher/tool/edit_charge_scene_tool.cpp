#include "schrodinger/sketcher/tool/edit_charge_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include <algorithm>

namespace schrodinger
{
namespace sketcher
{

EditChargeSceneTool::EditChargeSceneTool(ChargeTool charge_tool, Scene* scene,
                                         MolModel* mol_model) :
    SceneToolWithPredictiveHighlighting(scene, mol_model),
    m_charge_tool(charge_tool)
{
    m_highlight_types = InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
}

void EditChargeSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const RDKit::Atom* atom = atom_item->getAtom();
        auto new_charge = atom->getFormalCharge() +
                          (m_charge_tool == ChargeTool::INCREASE ? 1 : -1);
        // cap the charge between MIN_CHARGE and MAX_CHARGE
        new_charge = std::clamp(new_charge, MIN_ATOM_CHARGE, MAX_ATOM_CHARGE);
        m_mol_model->setAtomCharge(atom, new_charge);
    }
}

} // namespace sketcher
} // namespace schrodinger
