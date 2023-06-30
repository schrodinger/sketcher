#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

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
    AbstractSceneTool(scene, mol_model)
{
}

void MoveRotateSceneTool::onDragStart(QGraphicsSceneMouseEvent* event)
{
    std::vector<const RDKit::Atom*> selected_atoms;

    for (auto& item : m_scene->selectedItems()) {
        if (auto atom_item = dynamic_cast<AtomItem*>(item)) {
            selected_atoms.push_back(atom_item->getAtom());
        }
    }

    if (m_rotation_item.isInsideHandle(event->scenePos())) {
        m_action = Action::ROTATE;
        setAtomsToMove(selected_atoms);
    } else {
        m_action = Action::TRANSLATE;
        if (m_move_selection_item.rect().contains(event->scenePos())) {
            setAtomsToMove(selected_atoms);
        }
    }
}

void MoveRotateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    switch (m_action) {
        case Action::ROTATE:
            rotate(event, m_rotation_item.getPivotPoint(), m_atoms_to_move);
            break;
        case Action::TRANSLATE:
            translate(event, m_atoms_to_move);
            break;
        default:
            break;
    }
    updateMoveSelectionItem();
}

void MoveRotateSceneTool::onDragRelease(QGraphicsSceneMouseEvent* event)
{
    m_action = Action::NONE;
    setAtomsToMove({});
}

void MoveRotateSceneTool::setAtomsToMove(
    const std::vector<const RDKit::Atom*>& atoms)
{
    m_atoms_to_move = atoms;
}

void MoveRotateSceneTool::rotate(
    QGraphicsSceneMouseEvent* const event, QPointF pivot_point,
    const std::vector<const RDKit::Atom*>& atoms_to_move)
{
    // calculate angle increment of the event,
    // considering the centroid of the molecule as the origin
    QLineF current_position_line(pivot_point, event->scenePos());
    QLineF last_position_line(pivot_point, event->lastScenePos());
    auto angle = last_position_line.angleTo(current_position_line);
    m_mol_model->rotateByAngle(angle, to_mol_xy(pivot_point), atoms_to_move);
    rotateRotationItem(event);
}

void MoveRotateSceneTool::rotateRotationItem(QGraphicsSceneMouseEvent* event)
{
    auto center = m_rotation_item.getPivotPoint();
    auto prev_angle = QLineF(center, event->lastScenePos()).angle();
    auto angle = QLineF(center, event->scenePos()).angle();
    m_rotation_item.setArmAngle(m_rotation_item.getArmAngle() + angle -
                                prev_angle);
}

void MoveRotateSceneTool::translate(
    QGraphicsSceneMouseEvent* const event,
    const std::vector<const RDKit::Atom*>& atoms_to_move)
{
    auto distance = event->scenePos() - event->lastScenePos();
    m_mol_model->translateByVector(to_mol_xy(distance), atoms_to_move);
    m_rotation_item.setPivotPoint(m_rotation_item.getPivotPoint() + distance);
}

void MoveRotateSceneTool::updateRotationItem()
{
    if (m_mol_model->getMol()->getNumAtoms() == 0) {
        return;
    }
    m_rotation_item.setArmAngle(0);
    auto pivot_point = findPivotPointForRotation();
    m_rotation_item.setPivotPoint(to_scene_xy(pivot_point));

    // Show the rotation handle if there is more than one selected atom, or
    // if there is no selection and the molecule has more than one atom
    int num_atoms = (m_mol_model->getSelectedAtoms().size() > 0
                         ? m_mol_model->getSelectedAtoms().size()
                         : m_mol_model->getMol()->getNumAtoms());
    m_rotation_item.setDrawRotationHandle(num_atoms > 1);
}

void MoveRotateSceneTool::updateMoveSelectionItem()
{
    m_move_selection_item.setRect(m_scene->getSelectionRect());
}

std::vector<QGraphicsItem*> MoveRotateSceneTool::getGraphicsItems()
{
    updateRotationItem();
    updateMoveSelectionItem();

    return {&m_rotation_item, &m_move_selection_item};
}

RDGeom::Point3D MoveRotateSceneTool::findPivotPointForRotation()
{
    // if there's no coordinates, return the origin
    auto mol = m_mol_model->getMol();
    if (mol->getNumAtoms() == 0) {
        return RDGeom::Point3D(0, 0, 0);
    }

    auto selected_atoms = m_mol_model->getSelectedAtoms();
    if (selected_atoms.empty()) {
        return find_centroid(*mol);
    }

    // Find all bonds with one selected atom and one unselected atom
    std::unordered_set<const RDKit::Bond*> crossing_bonds;
    for (auto bond : mol->bonds()) {
        if (selected_atoms.count(bond->getBeginAtom()) !=
            selected_atoms.count(bond->getEndAtom())) {
            crossing_bonds.insert(bond);
        }
    }
    // If there is only one, use the selected atom as the pivot point
    if (crossing_bonds.size() == 1) {
        auto bond = *crossing_bonds.begin();
        auto begin_selected = selected_atoms.count(bond->getBeginAtom()) == 1;
        auto pivot_point_idx = begin_selected ? bond->getBeginAtom()->getIdx()
                                              : bond->getEndAtom()->getIdx();
        return mol->getConformer().getAtomPos(pivot_point_idx);
    }
    return find_centroid(*mol, selected_atoms);
}

} // namespace sketcher
} // namespace schrodinger