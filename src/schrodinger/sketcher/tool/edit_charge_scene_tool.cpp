#include "schrodinger/sketcher/tool/edit_charge_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include <algorithm>

namespace schrodinger
{
namespace sketcher
{

EditChargeSceneTool::EditChargeSceneTool(ChargeTool charge_tool, Scene* scene,
                                         MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_charge_tool(charge_tool)
{
    m_highlight_types = InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
}

void EditChargeSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const RDKit::Atom* atom = atom_item->getAtom();
        auto new_charge = atom->getFormalCharge() +
                          (m_charge_tool == ChargeTool::INCREASE ? 1 : -1);
        // cap the charge between MIN_CHARGE and MAX_CHARGE
        new_charge =
            std::clamp(new_charge, -ATOM_CHARGE_LIMIT, ATOM_CHARGE_LIMIT);
        m_mol_model->setAtomCharge(atom, new_charge);
    }
    StandardSceneToolBase::onLeftButtonClick(event);
}

QPixmap EditChargeSceneTool::createDefaultCursorPixmap() const
{
    QString path = m_charge_tool == ChargeTool::INCREASE
                       ? ":/icons/atom_charge_plus.svg"
                       : ":/icons/atom_charge_minus.svg";
    return cursor_hint_from_svg(path);
}

} // namespace sketcher
} // namespace schrodinger
