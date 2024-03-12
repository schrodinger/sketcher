#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"

#include <algorithm>

#include <GraphMol/MolOps.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger
{
namespace sketcher
{

MoveSelectionItem::MoveSelectionItem() : QGraphicsRectItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));
}

MoveRotateSceneTool::MoveRotateSceneTool(Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model)
{
    m_allow_context_menu = false;
    m_highlight_types = InteractiveItemFlag::ALL;
    updateRotationItem();
    updateMoveSelectionItem();
}

void MoveRotateSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    auto selected_atoms = m_mol_model->getSelectedAtoms();
    auto selected_non_mol_objs = m_mol_model->getSelectedNonMolecularObjects();

    QString description;
    if (m_rotation_item.isInsideHandle(event->scenePos())) {
        emit newCursorHintRequested(m_rotate_cursor_hint);
        m_action = Action::ROTATE;
        description = "Rotate";
        setObjectsToMove(selected_atoms, selected_non_mol_objs);
    } else {
        emit newCursorHintRequested(m_translate_cursor_hint);
        m_action = Action::TRANSLATE;
        description = "Translate";
        if (m_move_selection_item.rect().contains(event->scenePos())) {
            setObjectsToMove(selected_atoms, selected_non_mol_objs);
        }
    }
    // begin an undo macro in case we have to merge atoms at the end of the drag
    m_mol_model->beginUndoMacro(description);
    StandardSceneToolBase::onLeftButtonDragStart(event);
}

void MoveRotateSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    switch (m_action) {
        case Action::ROTATE:
            rotate(event, m_rotation_item.getPivotPoint(), m_atoms_to_move,
                   m_non_mol_objs_to_move);
            break;
        case Action::TRANSLATE:
            translate(event, m_atoms_to_move, m_non_mol_objs_to_move);
            break;
        default:
            break;
    }

    updateMergeHintItem();
    updateMoveSelectionItem();
    StandardSceneToolBase::onLeftButtonDragStart(event);
}

void MoveRotateSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    finishDrag();
    m_action = Action::NONE;
    m_rotation_item.setAngleValue("");
    updateMergeHintItem();
    StandardSceneToolBase::onLeftButtonDragMove(event);
}

void MoveRotateSceneTool::updateGraphicsItemsDuringRotation(
    QGraphicsSceneMouseEvent* event)
{
    auto center = m_rotation_item.getPivotPoint();
    auto prev_angle = QLineF(center, event->lastScenePos()).angle();
    auto angle = QLineF(center, event->scenePos()).angle();
    m_rotation_item.setArmAngle(m_rotation_item.getArmAngle() + angle -
                                prev_angle);

    auto angle_to_display =
        QLineF(center, event->scenePos())
            .angleTo(QLineF(center, m_drag_start_scene_pos));
    m_rotation_item.setAngleValue(QString::number(angle_to_display, 'f', 0));
}

void MoveRotateSceneTool::updateGraphicsItemsDuringTranslation(
    const QPointF& distance)
{
    m_rotation_item.setPivotPoint(m_rotation_item.getPivotPoint() + distance);
}

void MoveRotateSceneTool::updateRotationItem()
{
    if (m_mol_model->isEmpty()) {
        return;
    }
    m_rotation_item.setArmAngle(0);
    auto pivot_point = findPivotPointForRotation();
    m_rotation_item.setPivotPoint(to_scene_xy(pivot_point));

    // Show the rotation handle if there is more than one selected atom or
    // non-molecular object, or if there is no selection and the molecule has
    // more than one atom or non-molecular object
    auto sel_atoms = m_mol_model->getSelectedAtoms();
    auto sel_nmo = m_mol_model->getSelectedNonMolecularObjects();
    int num_atoms = sel_atoms.empty() && sel_nmo.empty()
                        ? m_mol_model->getMol()->getNumAtoms() +
                              m_mol_model->getNonMolecularObjects().size()
                        : sel_atoms.size() + sel_nmo.size();
    m_rotation_item.setDrawRotationHandle(num_atoms > 1);
}

void MoveRotateSceneTool::updateMoveSelectionItem()
{
    m_move_selection_item.setRect(m_scene->getSelectionRect());
}

std::vector<QGraphicsItem*> MoveRotateSceneTool::getGraphicsItems()
{
    auto graphics_items = StandardSceneToolBase::getGraphicsItems();
    graphics_items.push_back(&m_rotation_item);
    graphics_items.push_back(&m_move_selection_item);
    return graphics_items;
}

void MoveRotateSceneTool::onStructureUpdated()
{
    /**
     * don't update the rotation item if we're in the middle of
     * a drag othewise it will get reset to the original position.
     */
    if (!m_mouse_pressed) {
        updateRotationItem();
    }
    updateMoveSelectionItem();
}

QPixmap MoveRotateSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_move_rotate.svg");
}

} // namespace sketcher
} // namespace schrodinger